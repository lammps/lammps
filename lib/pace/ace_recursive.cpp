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

#include "ace_recursive.h"

#include "ace_abstract_basis.h"
#include "ace_types.h"

/* ------------------------------------------------------------
 *        ACEDAG Implementation
 *        (just the DAG construction, not the traversal)
 * ------------------------------------------------------------ */

/* Notes on different Tags: 
 *   rec1 - first basic implementation
 *   rec2 - avoid index arithmetic, contiguous layout, 
 *          canonical evaluator in ACE.jl format
 *   rec3 - split nodes into interior and leaf nodes
 */

void ACEDAG::init(Array2D<int> xAspec, 
                  Array2D<int> AAspec, 
                  Array1D<int> orders, 
                  Array2D<DOUBLE_TYPE> jl_coeffs,
                  int _heuristic ) {

    // remember which heuristic we want to use!
    heuristic = _heuristic;

    /* stage one of the graph is just extracting the A basis from 
     * the tensor product format into the linear list; all information 
     * for that is already stored in Aspec, and the only thing to do here is 
     * to construct zero-coefficients. Still we have to copy Aspec, since 
     * the one we have here will (may?) be deleted. */
    int num1 = xAspec.get_dim(0);
    Aspec = xAspec; //YL: just copying the multiarray: Aspec = xAspec;

    /* fill the one-particle basis into the DAGmap 
     *   DAGmap[ (i1,...,in) ] = iAA index where the (i1,...,in) basis functions 
     *   lives.
     */
    TDAGMAP DAGmap;
    for (int iA = 0; iA < num1; iA++) {
        vector<int> a(1);
        a[0] = iA;
        DAGmap[a] = iA; 
    }

    /* For stage 2 we now want to construct the actual recursion; the 
       recursion info will be stored in DAGspec, while the 
       coefficients go into DAGcoeffs. These arrays are initialized to 
       length `num2`, but they will have to grow as we add additional 
        artificial nodes into the graph. 
        
        initially we treat all nodes as having children, but then in a 
        second stage below we reorganize. */
    int num2 = AAspec.get_dim(0);
    int ndensity = jl_coeffs.get_dim(1);
    nodes_pre.resize(2*num2, 2);
    coeffs_pre.resize(2*num2, ndensity);
    
    /* the first basis function we construct will get index num1, 
     * since there are already num1 one-particle basis functions 
     * to collect during stage 1 */
    dag_idx = num1;
    /* main loop over AA basis set to transform into DAG  */
    for (int iAA = 0; iAA < num2; iAA++) {
        // create a vector representing the current basis function 
        int ord = orders(iAA);
        vector<int> aa(ord);
        for (int t = 0; t < ord; t++) aa[t] = AAspec(iAA, t); 
        vector<DOUBLE_TYPE> c(ndensity);
        for (int p = 0; p < ndensity; p++) c[p] = jl_coeffs(iAA, p);
        insert_node(DAGmap, aa, c);
    }

    /* convert to 3-stage format through reordering
     * interior nodes first, then leaf nodes  */

    // allocate storage 
    num_nodes = dag_idx;   // store total size of dag
    // num_nodes - num1 = number of many-body nodes.
    nodes.resize(num_nodes - num1, 2);
    coeffs.resize(num_nodes - num1, ndensity);

    // find out which nodes have children 
    haschild.resize(num_nodes - num1);
    haschild.fill(false);
    for (int iAA = 0; iAA < num_nodes - num1; iAA++) {
        if (nodes_pre(iAA, 0) >= num1)
            haschild(nodes_pre(iAA, 0)-num1) = true; 
        if (nodes_pre(iAA, 1) >= num1)
        haschild(nodes_pre(iAA, 1)-num1) = true; 
    }

    // to reorder the graph we need a fresh map from preordered indices  to 
    // postordered indices; for the 1-particle basis the order remains the same. 
    // TODO: doesn't need to be a map, could be a vector.
    map<int, int> neworder;
    for (int iA = 0; iA < num1; iA++) 
        neworder[iA] = iA;

    // insert all interior nodes 
    num2_int = 0;
    num2_leaf = 0;
    dag_idx = num1;
    int i1, i2, i1pre, i2pre; 
    for (int iAA = 0; iAA < num_nodes - num1; iAA++) {
        if (haschild(iAA)) {
            num2_int += 1;
            // indices into AAbuf before reordering
            i1pre = nodes_pre(iAA, 0);
            i2pre = nodes_pre(iAA, 1);
            // indices into AAbuf after reordering 
            i1 = neworder[i1pre];
            i2 = neworder[i2pre];
            // insert the current node : iAA is old order, dag_idx is new order
            neworder[num1+iAA] = dag_idx; 
            nodes(dag_idx-num1, 0) = i1; 
            nodes(dag_idx-num1, 1) = i2; 
            for (int t = 0; t < ndensity; t++)
                coeffs(dag_idx-num1, t) = coeffs_pre(iAA, t);            
            dag_idx++;
        }
    }

    // insert all leaf nodes
    for (int iAA = 0; iAA < num_nodes - num1; iAA++) {
        if (!haschild(iAA)) {
            num2_leaf += 1;
            // indices into AAbuf before reordering
            i1pre = nodes_pre(iAA, 0);
            i2pre = nodes_pre(iAA, 1);            
            // insert the current node : no need to remember the new order now
            nodes(dag_idx-num1, 0) = neworder[i1pre];
            nodes(dag_idx-num1, 1) = neworder[i2pre];
            for (int t = 0; t < ndensity; t++)
                coeffs(dag_idx-num1, t) = coeffs_pre(iAA, t);
            dag_idx++;
        }
    }
#ifdef DEBUG
    cout << "num2_int = " << num2_int << "; num2_leaf = " << num2_leaf << "\n";
#endif
    // free up memory that is no longer needed
    nodes_pre.resize(0, 0);
    coeffs_pre.resize(0, 0);
    haschild.resize(0);

    /* finalize dag: allocate buffer storage  */
    AAbuf.resize(num1 + num2_int);
    w.resize(num_nodes);   
    // TODO: technically only need num1 + num2_int for w, this can save  one
    //       memory access later, probably not worth the crazy code duplication.
}

void ACEDAG::insert_node(TDAGMAP &DAGmap, vector<int> a, vector<DOUBLE_TYPE> c) {
    /* start with a list of all possible partitions into 2 groups 
     * and check whether any of these nodes are already in the dag */ 
    auto partitions = find_2partitions(a);
    int ndensity = c.size();
    int num1 = get_num1();

    // TODO: first try to find partitions into nodes that are already parents
    //       that way we will get more leaf nodes!
    for (TPARTITION const& p : partitions) {
        /* this is the good case; the parent nodes are both already in the 
         * graph; add the new node and return. This is also the only place in the 
         * code where an actual insert happens. */
        if (DAGmap.count(p.first) && DAGmap.count(p.second)) {
            if (nodes_pre.get_dim(0) < dag_idx + 1) { //check if array is sufficiently large
                int newsize = (dag_idx * 3) / 2;
                nodes_pre.resize(newsize, 2); // grow arrays if necessary
                coeffs_pre.resize(newsize, ndensity);
            }
            int i1 = DAGmap[p.first]; 
            int i2 = DAGmap[p.second];
            nodes_pre(dag_idx - num1, 0) = i1;
            nodes_pre(dag_idx - num1, 1) = i2;
            DAGmap[a] = dag_idx; 
            for (int p = 0; p < ndensity; p++)  
                coeffs_pre(dag_idx - num1, p) = c[p];
            dag_idx += 1;
            return;
        }
    }

    /* if we are here, then this means, the new node cannot yet be inserted. 
     * We first need to insert some intermediate auxiliary nodes.  For this 
     * we use a simple heuristic: 
     *    try to find a partition where one of the two nodes are already 
     *    in the graph, if there are several, then we remember the longest
     *    (this is a very greedy heuristic!!)
     *  .... (continue below) ....
     */
    TPARTITION longest;
    int longest_length = 0; 
    for (auto const& p : partitions) {
        int len = 0;
        if (DAGmap.count(p.first)) {
            len = p.first.size();
        } else if (DAGmap.count(p.second)) {
            len = p.second.size();
        }
        if ((len > 0) && (len > longest_length)) {
            longest_length = len;
            longest = p;
        }
    }

    /* sanity check */
    if (longest_length == 0) {
        std::stringstream error_message;
        error_message << "WARNING : something has gone horribly wrong! `longest_length == 0`! \n";
        error_message << "a = [";
        for (int t = 0; t < a.size(); t++)
            error_message << a[t] << ", ";
        error_message << "]\n";
        throw std::logic_error(error_message.str());
//        return;
    }

    /* If there is a partition with one component already in the graph, 
     * then we only need to add in the other component. Note that there will 
     * always be at least one such partition, namely all those containing 
     * a one-element node e.g. (1,2,3,4) -> (1,) (2,3,4)  then (1,) is 
     * a one-particle basis function and hence always in the graph. 
     * If heuristic == 0, then we just take one of those partitionas and move on.
     * 
     * We also accept the found partition if longest_length > 1.
     * And we also accept it if we have a 2- or 3-correlation. 
     */

    if (     (heuristic == 0) 
          || (longest_length > 1)
          || (a.size() <= 3))       {
        /* insert the other node that isn't in the DAG yet 
        * this is an artificial node so it gets zero-coefficients 
        * This step is recursive, so more than one node might be inserted here */
        vector<DOUBLE_TYPE> cz(ndensity); 
        for (int i = 0; i < ndensity; i++) cz[i] = 0.0; 
        TPARTITION p = longest;
        if (DAGmap.count(p.first)) 
            insert_node(DAGmap, p.second, cz);
        else
            insert_node(DAGmap, p.first, cz);
    }

    /* Second heuristic : heuristic == 1
     * Focus on inserting artificial 2-correlations
     */
    else if (heuristic == 1) {
        // and we also know that longest_length == 1 and nu = a.size >= 4.
        int nu = a.size();
        // generate an artificial partition
        vector<int> a1(2);
        for (int i = 0; i < 2; i++) a1[i] = a[i];
        vector<int> a2(nu - 2);
        for (int i = 0; i < nu - 2; i++) a2[i] = a[2 + i];
        vector<DOUBLE_TYPE> cz(ndensity);
        for (int i = 0; i < cz.size(); i++) cz[i] = 0.0;
        // and insert both (we know neither are in the DAG yet)
        insert_node(DAGmap, a1, cz);
        insert_node(DAGmap, a2, cz);
    } else {
        cout << "WARNING : something has gone horribly wrong! \n";
        //  TODO: Throw and error here?!?
        return;
    }
    


    /* now we should be ready to insert the entire tuple `a` since there is now 
     * an eligible parent pair. Here we recompute the partition of `a`, but 
     * that's a small price to pay for a clearer code. Maybe this can be 
     * optimized a bit by wrapping it all into a while loop or having a second
     * version of `insert_node` ... */
    insert_node(DAGmap, a, c);
}

TPARTITIONS ACEDAG::find_2partitions(vector<int> v) {
    int N = v.size(); 
    int zo;
    TPARTITIONS partitions;
    TPARTITION part;
    /* This is a fun little hack to extract all subsets of the indices 1:N 
     * the number i will have binary representation with each digit indicating 
     * whether or not that index belongs to the selected subset */
    for (int i = 1; i < (1<<N)/2; i++){
        int N1 = 0, N2 = 0; 
        int p = 1; 
        for (int n = 0; n < N; n++) {
            zo = ((i / p) % 2);
            N1 += zo; 
            N2 += 1-zo;
            p *= 2; 
        }
        /* convert to a more useful representation in terms of vector */
        vector<int> v1(N1);
        vector<int> v2(N2);
        int i1 =0, i2 = 0;
        p = 1;
        for (int n = 0; n < N; n++) {
            zo = ((i / p) % 2);
            p *= 2;
            if (zo == 1) { 
                v1[i1] = v[n]; 
                i1 += 1;
            } else {
                v2[i2] = v[n]; 
                i2 += 1;
            }
        }
        part = make_pair(v1, v2);
        partitions.push_back(part);
    }
    return partitions;
}

void ACEDAG::print() {
    cout << "DAG Specification: \n" ;
    cout << "          n1 : " << get_num1() << "\n"; 
    cout << "          n2 : " << get_num2() << "\n"; 
    cout << "   num_nodes : " << num_nodes << "\n"; 
    cout << "--------------------\n";
    cout << "A-spec: \n";
    for (int iA = 0; iA < get_num1(); iA++) {
        cout << iA << " : " << Aspec(iA, 0) << 
            Aspec(iA, 1) << Aspec(iA, 2) << Aspec(iA, 3) << "\n";
    }

    cout << "-----------\n";
    cout << "AA-tree\n";
    
    for (int iAA = 0; iAA < get_num2(); iAA++) {
        cout << iAA + get_num1() << " : " << 
            nodes(iAA, 0) << ", " << nodes(iAA, 1) << "\n";
    }
}


/* ------------------------------------------------------------
 *        ACERecursiveEvaluator
 * ------------------------------------------------------------ */


void ACERecursiveEvaluator::set_basis(ACECTildeBasisSet &bas, int heuristic) {
    basis_set = &bas;
    init(basis_set, heuristic);
}

void ACERecursiveEvaluator::init(ACECTildeBasisSet *basis_set, int heuristic) {

    ACEEvaluator::init(basis_set);


    weights.init(basis_set->nelements, basis_set->nradmax + 1, basis_set->lmax + 1,
                 "weights");

    weights_rank1.init(basis_set->nelements, basis_set->nradbase, "weights_rank1");


    DG_cache.init(1, basis_set->nradbase, "DG_cache");
    DG_cache.fill(0);

    R_cache.init(1, basis_set->nradmax, basis_set->lmax + 1, "R_cache");
    R_cache.fill(0);

    DR_cache.init(1, basis_set->nradmax, basis_set->lmax + 1, "DR_cache");
    DR_cache.fill(0);

    Y_cache.init(1, basis_set->lmax + 1, "Y_cache");
    Y_cache.fill({0, 0});

    DY_cache.init(1, basis_set->lmax + 1, "dY_dense_cache");
    DY_cache.fill({0.});

    //hard-core repulsion
    DCR_cache.init(1, "DCR_cache");
    DCR_cache.fill(0);
    dB_flatten.init(basis_set->max_dB_array_size, "dB_flatten");

    /* convert to ACE.jl format to prepare for construction of DAG
     * This will fill the arrays jl_Aspec, jl_AAspec, jl_orders
     */
    acejlformat();

    // test_acejlformat(); 
    
    // now pass this info into the DAG 
    dag.init(jl_Aspec, jl_AAspec, jl_orders, jl_coeffs, heuristic);

    // finally empty the temporary arrays to clear up the memory...
    // TODO 
}


void ACERecursiveEvaluator::acejlformat() {

    int func_ms_ind = 0;
    int func_ms_t_ind = 0;// index for dB
    int j, jj, func_ind, ms_ind;

    SPECIES_TYPE mu_i = 0;//TODO: multispecies
    const SHORT_INT_TYPE total_basis_size = basis_set->total_basis_size[mu_i];
    ACECTildeBasisFunction *basis = basis_set->basis[mu_i];

    int AAidx = 0;
    RANK_TYPE order, t;
    SPECIES_TYPE *mus;
    NS_TYPE *ns;
    LS_TYPE *ls;
    MS_TYPE *ms;

    /* transform basis into new format: 
       [A1 ... A_num1]
       [(i1,i2)(i1,i2)(...)(i1,i2,i3)(...)]
       where each ia represents an A_{ia}
    */

    /* compute max values for mu, n, l, m */
    SPECIES_TYPE maxmu = 0; //TODO: multispecies
    NS_TYPE maxn = basis_set->nradmax;
    LS_TYPE maxl = basis_set->lmax;
    RANK_TYPE maxorder = basis_set->rankmax;
    const DENSITY_TYPE ndensity = basis_set->ndensitymax;

    int num1 = 0;


    /* create a 4D lookup table for the 1-p basis 
     * TODO: replace with a map??
     */
    Array4D<int> A_lookup(int(maxmu+1), int(maxn), int(maxl+1), int(2*maxl+1));
    for (int mu = 0; mu < maxmu+1; mu++) 
        for (int n = 0; n < maxn; n++)
            for (int l = 0; l < maxl+1; l++)
                for (int m = 0; m < 2*maxl+1; m++)
                    A_lookup(mu, n, l, m) = -1;
    int A_idx = 0;  // linear index of A basis function (1-particle)
    for (func_ind = 0; func_ind < total_basis_size; ++func_ind) {
        ACECTildeBasisFunction *func = &basis[func_ind];
//        func->print();
        order = func->rank; mus = func->mus; ns = func->ns; ls = func->ls;
        for (ms_ind = 0; ms_ind < func->num_ms_combs; ++ms_ind, ++func_ms_ind) {
            ms = &func->ms_combs[ms_ind * order]; 
            for (t = 0; t < order; t++) {
                int iA = A_lookup(mus[t], ns[t]-1, ls[t], ms[t]+ls[t]);
                if (iA == -1) {
                    A_lookup(mus[t], ns[t] - 1, ls[t], ms[t] + ls[t]) = A_idx;
                    A_idx += 1;
                }
            }
        }
    }

    /* create the reverse list: linear indixes to mu,l,m,n
       this keeps only the basis functions we really need */
    num1 = A_idx;
    Array2D<int> & Aspec = jl_Aspec; 
    Aspec.resize(num1, 4);
    // Array2D<int> Aspec(num1, 4);
    for (int mu = 0; mu <= maxmu; mu++)
        for (int n = 1; n <= maxn; n++)
            for (int l = 0; l <= maxl; l++)
                for (int m = -l; m <= l; m++) {
                    int iA = A_lookup(mu, n-1, l, l+m);
                    if (iA != -1) {
                        Aspec(iA, 0) = mu;
                        Aspec(iA, 1) = n;
                        Aspec(iA, 2) = l;
                        Aspec(iA, 3) = m;
                    }
                }

    /* ============ HALF-BASIS TRICK START ============ */
    for (func_ind = 0; func_ind < total_basis_size; ++func_ind) {
        ACECTildeBasisFunction *func = &basis[func_ind];
        order = func->rank; mus = func->mus; ns = func->ns; ls = func->ls;
        if (!( (mus[0] <= maxmu) && (ns[0] <= maxn) && (ls[0] <= maxl) ))
            continue; 

        for (ms_ind = 0; ms_ind < func->num_ms_combs; ++ms_ind, ++func_ms_ind) {
            ms = &func->ms_combs[ms_ind * order]; 

            // find first positive and negative index
            int pos_idx = order + 1;
            int neg_idx = order + 1;
            for (t = order-1; t >= 0; t--)
                if (ms[t] > 0) pos_idx = t; 
                else if (ms[t] < 0) neg_idx = t; 

            // if neg_idx < pos_idx then this means that ms is non-zero 
            // and that the first non-zero index is negative, hence this is 
            // a negative-sign tuple which we want to combine into 
            // its opposite.
            if (neg_idx < pos_idx) {
                // find the opposite tuple
                int ms_ind2 = 0;
                MS_TYPE *ms2;
                bool found_opposite = false; 
                for (ms_ind2 = 0; ms_ind2 < func->num_ms_combs; ++ms_ind2) {
                    ms2 = &func->ms_combs[ms_ind2 * order]; 
                    bool isopposite = true; 
                    for (t = 0; t < order; t++)
                        if (ms[t] != -ms2[t]) {
                            isopposite = false; 
                            break;
                        }
                    if (isopposite) {
                        found_opposite = true;
                        break;
                    }
                }

                if (ms_ind == ms_ind2) {
                    cout << "WARNING - ms_ind == ms_ind2 \n";
                }

                // now we need to overwrite the coefficients
                if (found_opposite)  {
                    int sig = 1;
                    for (t = 0; t < order; t++)
                        if (ms[t] < 0)
                            sig *= -1;
                    for (int p = 0; p < ndensity; ++p) {
                        func->ctildes[ms_ind2 * ndensity + p] += 
                                func->ctildes[ms_ind * ndensity + p];
                        func->ctildes[ms_ind * ndensity + p] = 0.0;
                    }
                }
            }
        }
    }

    // /* ============ HALF-BASIS TRICK END ============ */


    /* count number of basis functions, keep only non-zero!!  */
    int num2 = 0;
    for (func_ind = 0; func_ind < total_basis_size; ++func_ind)  {
        ACECTildeBasisFunction *func = &basis[func_ind];
        for (ms_ind = 0; ms_ind < (&basis[func_ind])->num_ms_combs; ++ms_ind, ++func_ms_ind) {
            // check that the coefficients are actually non-zero 
            bool isnonzero = false; 
            for (DENSITY_TYPE p = 0; p < ndensity; ++p)
                if (func->ctildes[ms_ind * ndensity + p] != 0.0)
                    isnonzero = true;
            if (isnonzero)
                num2++;
        }
    }


    /* Now create the AA basis links into the A basis */
    num1 = A_idx;   // total number of A-basis functions that we keep
    // Array1D<int> AAorders(num2);
    Array1D<int> & AAorders = jl_orders; 
    AAorders.resize(num2);
    // Array2D<int> AAspec(num2, maxorder);    // specs of AA basis functions
    Array2D<int> & AAspec = jl_AAspec; 
    AAspec.resize(num2, maxorder);
    jl_coeffs.resize(num2, ndensity);
    AAidx = 0;                          // linear index into AA basis function 
    int len_flat = 0;
    for (func_ind = 0; func_ind < total_basis_size; ++func_ind) {
        ACECTildeBasisFunction *func = &basis[func_ind];
        order = func->rank; mus = func->mus; ns = func->ns; ls = func->ls;
        if (!((mus[0] <= maxmu) && (ns[0] <= maxn) && (ls[0] <= maxl)))        //fool-proofing of functions
            continue; 

        for (ms_ind = 0; ms_ind < func->num_ms_combs; ++ms_ind, ++func_ms_ind) {
            ms = &func->ms_combs[ms_ind * order]; 

            // check that the coefficients are actually non-zero 
            bool iszero = true; 
            for (DENSITY_TYPE p = 0; p < ndensity; ++p)
                if (func->ctildes[ms_ind * ndensity + p] != 0.0)
                    iszero = false;
            if (iszero) continue;

            AAorders(AAidx) = order;
            for (t = 0; t < order; t++) {
                int Ait = A_lookup(int(mus[t]), int(ns[t]-1), int(ls[t]), int(ms[t])+int(ls[t]));
                AAspec(AAidx, t) = Ait; 
                len_flat += 1;
            }
            for (t = order; t < maxorder; t++) AAspec(AAidx, t) = -1; 
            /* copy over the coefficients */
            for (DENSITY_TYPE p = 0; p < ndensity; ++p)
                jl_coeffs(AAidx, p) = func->ctildes[ms_ind * ndensity + p];
            AAidx += 1;
        }
    }

    // flatten the AAspec array 
    jl_AAspec_flat.resize(len_flat);
    int idx_spec = 0;
    for (int AAidx = 0; AAidx < jl_AAspec.get_dim(0); AAidx++) 
        for (int p = 0; p < jl_orders(AAidx); p++, idx_spec++)
            jl_AAspec_flat(idx_spec) = jl_AAspec(AAidx, p);

}

void ACERecursiveEvaluator::test_acejlformat() {

    Array2D<int> AAspec = jl_AAspec;
    Array2D<int> Aspec = jl_Aspec;
    Array1D<int> AAorders = jl_orders;
    cout << "num2 = " << AAorders.get_dim(0) << "\n";
    int func_ms_ind = 0;
    int func_ms_t_ind = 0;// index for dB
    int j, jj, func_ind, ms_ind;

    SPECIES_TYPE mu_i = 0;
    const SHORT_INT_TYPE total_basis_size = basis_set->total_basis_size[mu_i];
    ACECTildeBasisFunction *basis = basis_set->basis[mu_i];

    RANK_TYPE order, t;
    SPECIES_TYPE *mus;
    NS_TYPE *ns;
    LS_TYPE *ls;
    MS_TYPE *ms;

    /* ==== test by printing the basis spec ====*/
    // TODO: convert this into an automatic consistency test 
    int iAA = 0;
    for (func_ind = 0; func_ind < total_basis_size; ++func_ind) {
        ACECTildeBasisFunction *func = &basis[func_ind];
        order = func->rank; mus = func->mus; ns = func->ns; ls = func->ls;
        // func->print();
        //loop over {ms} combinations in sum
        for (ms_ind = 0; ms_ind < func->num_ms_combs; ++ms_ind, ++func_ms_ind) {
            ms = &func->ms_combs[ms_ind * order]; 

            
            cout << iAA << " : |";
            for (t = 0; t < order; t++)
                cout << mus[t] << ";" << ns[t] << "," << ls[t] << "," << ms[t] << "|";
            cout << "\n"; 

            cout << "      ["; 
            for (t = 0; t < AAorders(iAA); t++) 
                cout << AAspec(iAA, int(t)) << ",";
            cout << "]\n";
            cout << "      |"; 
            for (t = 0; t < AAorders(iAA); t++)  { 
                int iA = AAspec(iAA, t);
                // cout << iA << ",";
                cout << Aspec(iA, 0) << ";" 
                     << Aspec(iA, 1) << "," 
                     << Aspec(iA, 2) << "," 
                     << Aspec(iA, 3) << "|";
            }
            cout << "\n";
            iAA += 1;
        }
    }
    /* ==== END TEST ==== */


}




void ACERecursiveEvaluator::resize_neighbours_cache(int max_jnum) {
    if(basis_set== nullptr) {
        throw std::invalid_argument("ACERecursiveEvaluator: basis set is not assigned");
    }
    if (R_cache.get_dim(0) < max_jnum) {

        //TODO: implement grow
        R_cache.resize(max_jnum, basis_set->nradmax, basis_set->lmax + 1);
        R_cache.fill(0);

        DR_cache.resize(max_jnum, basis_set->nradmax, basis_set->lmax + 1);
        DR_cache.fill(0);

        DG_cache.resize(max_jnum, basis_set->nradbase);
        DG_cache.fill(0);

        Y_cache.resize(max_jnum, basis_set->lmax + 1);
        Y_cache.fill({0, 0});

        DY_cache.resize(max_jnum, basis_set->lmax + 1);
        DY_cache.fill({0});

        //hard-core repulsion
        DCR_cache.init(max_jnum, "DCR_cache");
        DCR_cache.fill(0);
    }
}



// double** r - atomic coordinates of atom I
// int* types - atomic types if atom I
// int **firstneigh -  ptr to 1st J int value of each I atom. Usage: jlist = firstneigh[i];
// Usage: j = jlist_of_i[jj];
// jnum - number of J neighbors for each I atom.  jnum = numneigh[i];

void
ACERecursiveEvaluator::compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, const int jnum, const int *jlist) {
    if(basis_set== nullptr) {
        throw std::invalid_argument("ACERecursiveEvaluator: basis set is not assigned");
    }
    per_atom_calc_timer.start();
#ifdef PRINT_MAIN_STEPS
    printf("\n ATOM: ind = %d r_norm=(%f, %f, %f)\n",i, x[i][0], x[i][1], x[i][2]);
#endif
    DOUBLE_TYPE evdwl = 0, evdwl_cut = 0, rho_core = 0;
    DOUBLE_TYPE r_norm;
    DOUBLE_TYPE xn, yn, zn, r_xyz;
    DOUBLE_TYPE R, GR, DGR, R_over_r, DR, DCR;
    DOUBLE_TYPE *r_hat;

    SPECIES_TYPE mu_j;
    RANK_TYPE r, rank, t;
    NS_TYPE n;
    LS_TYPE l;
    MS_TYPE m, m_t;

    SPECIES_TYPE *mus;
    NS_TYPE *ns;
    LS_TYPE *ls;
    MS_TYPE *ms;

    int j, jj, func_ind, ms_ind;
    SHORT_INT_TYPE factor;

    ACEComplex Y{0}, Y_DR{0.};
    ACEComplex B{0.};
    ACEComplex dB{0};
    ACEComplex A_cache[basis_set->rankmax];
    
    ACEComplex dA[basis_set->rankmax];
    int spec[basis_set->rankmax];

    dB_flatten.fill({0.});

    ACEDYcomponent grad_phi_nlm{0}, DY{0.};

    //size is +1 of max to avoid out-of-boundary array access in double-triangular scheme
    ACEComplex A_forward_prod[basis_set->rankmax + 1];
    ACEComplex A_backward_prod[basis_set->rankmax + 1];

    DOUBLE_TYPE inv_r_norm;
    DOUBLE_TYPE r_norms[jnum];
    DOUBLE_TYPE inv_r_norms[jnum];
    DOUBLE_TYPE rhats[jnum][3];//normalized vector
    SPECIES_TYPE elements[jnum];
    const DOUBLE_TYPE xtmp = x[i][0];
    const DOUBLE_TYPE ytmp = x[i][1];
    const DOUBLE_TYPE ztmp = x[i][2];
    DOUBLE_TYPE f_ji[3];

    bool is_element_mapping = element_type_mapping.get_size() > 0;
    SPECIES_TYPE mu_i;
    if (is_element_mapping)
        mu_i = element_type_mapping(type[i]);
    else
        mu_i = type[i];

    const SHORT_INT_TYPE total_basis_size_rank1 = basis_set->total_basis_size_rank1[mu_i];
    const SHORT_INT_TYPE total_basis_size = basis_set->total_basis_size[mu_i];

    ACECTildeBasisFunction *basis_rank1 = basis_set->basis_rank1[mu_i];
    ACECTildeBasisFunction *basis = basis_set->basis[mu_i];

    DOUBLE_TYPE rho_cut, drho_cut, fcut, dfcut;
    DOUBLE_TYPE dF_drho_core;

    //TODO: lmax -> lmaxi (get per-species type)
    const LS_TYPE lmaxi = basis_set->lmax;

    //TODO: nradmax -> nradiali (get per-species type)
    const NS_TYPE nradiali = basis_set->nradmax;

    //TODO: nradbase -> nradbasei (get per-species type)
    const NS_TYPE nradbasei = basis_set->nradbase;

    //TODO: get per-species type number of densities
    const DENSITY_TYPE ndensity= basis_set->ndensitymax;

    neighbours_forces.resize(jnum, 3);
    neighbours_forces.fill(0);

    //TODO: shift nullifications to place where arrays are used
    weights.fill({0});
    weights_rank1.fill(0);
    A.fill({0});
    A_rank1.fill(0);
    rhos.fill(0);
    dF_drho.fill(0);

#ifdef EXTRA_C_PROJECTIONS
    basis_projections_rank1.init(total_basis_size_rank1, ndensity, "c_projections_rank1");
    basis_projections.init(total_basis_size, ndensity, "c_projections");
#endif

    //proxy references to spherical harmonics and radial functions arrays
    const Array2DLM<ACEComplex> &ylm = basis_set->spherical_harmonics.ylm;
    const Array2DLM<ACEDYcomponent> &dylm = basis_set->spherical_harmonics.dylm;

    const Array2D<DOUBLE_TYPE> &fr = basis_set->radial_functions->fr;
    const Array2D<DOUBLE_TYPE> &dfr = basis_set->radial_functions->dfr;

    const Array1D<DOUBLE_TYPE> &gr = basis_set->radial_functions->gr;
    const Array1D<DOUBLE_TYPE> &dgr = basis_set->radial_functions->dgr;

    loop_over_neighbour_timer.start();

    int jj_actual = 0;
    SPECIES_TYPE type_j = 0;
    int neighbour_index_mapping[jnum]; // jj_actual -> jj
    //loop over neighbours, compute distance, consider only atoms within with r<cutoff(mu_i, mu_j)
    for (jj = 0; jj < jnum; ++jj) {

        j = jlist[jj];
        xn = x[j][0] - xtmp;
        yn = x[j][1] - ytmp;
        zn = x[j][2] - ztmp;
        type_j = type[j];
        if (is_element_mapping)
            mu_j = element_type_mapping(type_j);
        else
            mu_j = type_j;

        DOUBLE_TYPE current_cutoff = basis_set->radial_functions->cut(mu_i, mu_j);
        r_xyz = sqrt(xn * xn + yn * yn + zn * zn);

        if (r_xyz >= current_cutoff)
            continue;

        inv_r_norm = 1 / r_xyz;

        r_norms[jj_actual] = r_xyz;
        inv_r_norms[jj_actual] = inv_r_norm;
        rhats[jj_actual][0] = xn * inv_r_norm;
        rhats[jj_actual][1] = yn * inv_r_norm;
        rhats[jj_actual][2] = zn * inv_r_norm;
        elements[jj_actual] = mu_j;
        neighbour_index_mapping[jj_actual] = jj;
        jj_actual++;
    }

    int jnum_actual = jj_actual;

    //ALGORITHM 1: Atomic base A
    for (jj = 0; jj < jnum_actual; ++jj) {
        r_norm = r_norms[jj];
        mu_j = elements[jj];
        r_hat = rhats[jj];

        //proxies
        Array2DLM<ACEComplex> &Y_jj = Y_cache(jj);
        Array2DLM<ACEDYcomponent> &DY_jj = DY_cache(jj);


        basis_set->radial_functions->evaluate(r_norm, basis_set->nradbase, nradiali, mu_i, mu_j);
        basis_set->spherical_harmonics.compute_ylm(r_hat[0], r_hat[1], r_hat[2], lmaxi);
        //loop for computing A's
        //rank = 1
        for (n = 0; n < basis_set->nradbase; n++) {
            GR = gr(n);
#ifdef DEBUG_ENERGY_CALCULATIONS
            printf("-neigh atom %d\n", jj);
            printf("gr(n=%d)(r=%f) = %f\n", n, r_norm, gr(n));
            printf("dgr(n=%d)(r=%f) = %f\n", n, r_norm, dgr(n));
#endif
            DG_cache(jj, n) = dgr(n);
            A_rank1(mu_j, n) += GR * Y00;
        }
        //loop for computing A's
        // for rank > 1
        for (n = 0; n < nradiali; n++) {
            auto &A_lm = A(mu_j, n);
            for (l = 0; l <= lmaxi; l++) {
                R = fr(n, l);
#ifdef DEBUG_ENERGY_CALCULATIONS
                printf("R(nl=%d,%d)(r=%f)=%f\n", n + 1, l, r_norm, R);
#endif

                DR_cache(jj, n, l) = dfr(n, l);
                R_cache(jj, n, l) = R;

                for (m = 0; m <= l; m++) {
                    Y = ylm(l, m);
#ifdef DEBUG_ENERGY_CALCULATIONS
                    printf("Y(lm=%d,%d)=(%f, %f)\n", l, m, Y.real, Y.img);
#endif
                    A_lm(l, m) += R * Y; //accumulation sum over neighbours
                    Y_jj(l, m) = Y;
                    DY_jj(l, m) = dylm(l, m);
                }
            }
        }

        //hard-core repulsion
        rho_core += basis_set->radial_functions->cr;
        DCR_cache(jj) = basis_set->radial_functions->dcr;

    } //end loop over neighbours

    //complex conjugate A's (for NEGATIVE (-m) terms)
    // for rank > 1
    for (mu_j = 0; mu_j < basis_set->nelements; mu_j++) {
        for (n = 0; n < nradiali; n++) {
            auto &A_lm = A(mu_j, n);
            for (l = 0; l <= lmaxi; l++) {
                //fill in -m part in the outer loop using the same m <-> -m symmetry as for Ylm
                for (m = 1; m <= l; m++) {
                    factor = m % 2 == 0 ? 1 : -1;
                    A_lm(l, -m) = A_lm(l, m).conjugated() * factor;
                }
            }
        }
    }    //now A's are constructed
    loop_over_neighbour_timer.stop();

    // ==================== ENERGY ====================

    energy_calc_timer.start();
#ifdef EXTRA_C_PROJECTIONS
    basis_projections_rank1.fill(0);
    basis_projections.fill(0);
#endif

    //ALGORITHM 2: Basis functions B with iterative product and density rho(p) calculation
    //rank=1
    for (int func_rank1_ind = 0; func_rank1_ind < total_basis_size_rank1; ++func_rank1_ind) {
        ACECTildeBasisFunction *func = &basis_rank1[func_rank1_ind];
//        ndensity = func->ndensity;
#ifdef PRINT_LOOPS_INDICES
        printf("Num density = %d r = 0\n",(int) ndensity );
        print_C_tilde_B_basis_function(*func);
#endif
        double A_cur = A_rank1(func->mus[0], func->ns[0] - 1);
#ifdef DEBUG_ENERGY_CALCULATIONS
        printf("A_r=1(x=%d, n=%d)=(%f)\n", func->mus[0], func->ns[0], A_cur);
        printf("     coeff[0] = %f\n", func->ctildes[0]);
#endif
        for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
            //for rank=1 (r=0) only 1 ms-combination exists (ms_ind=0), so index of func.ctildes is 0..ndensity-1
            rhos(p) += func->ctildes[p] * A_cur;
#ifdef EXTRA_C_PROJECTIONS
            //aggregate C-projections separately
            basis_projections_rank1(func_rank1_ind, p)+= func->ctildes[p] * A_cur;
#endif
        }
    } // end loop for rank=1

    // ================ START RECURSIVE EVALUATOR ====================
    // (rank > 1 only)

    /* STAGE 1: 
     * 1-particle basis is already evaluated, so we only need to 
     * copy it into the AA value buffer 
     */ 
    int num1 = dag.get_num1();
    for (int idx = 0; idx < num1; idx++) 
        dag.AAbuf(idx) = A( dag.Aspec(idx, 0), 
                            dag.Aspec(idx, 1)-1, 
                            dag.Aspec(idx, 2), 
                            dag.Aspec(idx, 3) );


    if (recursive) {
        /* STAGE 2: FORWARD PASS
        * Forward pass: go through the dag and store all intermediate results
        */

        // rhos.fill(0); note the rhos are already reset and started filling above!
        ACEComplex AAcur{0.0};
        int i1, i2;

        int * dag_nodes = dag.nodes.get_data();
        int idx_nodes = 0;

        DOUBLE_TYPE * dag_coefs = dag.coeffs.get_data();
        int idx_coefs = 0;
        
        int num2_int = dag.get_num2_int();
        int num2_leaf = dag.get_num2_leaf();

        // interior nodes (save AA)
        for (int idx = num1; idx < num1+num2_int; idx++) {    
            i1 = dag_nodes[idx_nodes]; idx_nodes++;
            i2 = dag_nodes[idx_nodes]; idx_nodes++;
            AAcur = dag.AAbuf(i1) * dag.AAbuf(i2);
            dag.AAbuf(idx) = AAcur;
            for (int p = 0; p < ndensity; p++, idx_coefs++)            
                rhos(p) += AAcur.real_part_product(dag_coefs[idx_coefs]);
        }

        // leaf nodes -> no need to store in AAbuf
        DOUBLE_TYPE AAcur_re = 0.0;
        for (int _idx = 0; _idx < num2_leaf; _idx++) {    
            i1 = dag_nodes[idx_nodes]; idx_nodes++;
            i2 = dag_nodes[idx_nodes]; idx_nodes++;
            AAcur_re = dag.AAbuf(i1).real_part_product(dag.AAbuf(i2));
            for (int p = 0; p < ndensity; p++, idx_coefs++)            
                rhos(p) += AAcur_re * dag_coefs[idx_coefs];
        }

    } else {

        /* non-recursive Julia-style evaluator implementation */
        // TODO: fix array access to enable bounds checking again???
        ACEComplex AAcur{1.0};
        int *AAspec = jl_AAspec_flat.get_data();
        DOUBLE_TYPE *coeffs = jl_coeffs.get_data();
        int idx_spec = 0;
        int idx_coefs = 0;
        int order = 0;
        int max_order = jl_AAspec.get_dim(1);
        for (int iAA = 0; iAA < jl_AAspec.get_dim(0); iAA ++) {
            AAcur = 1.0;
            order = jl_orders(iAA);
            for (int r = 0; r < order; r++, idx_spec++)
                AAcur *= dag.AAbuf( AAspec[idx_spec] );
            for (int p = 0; p < ndensity; p++, idx_coefs++)            
                rhos(p) += AAcur.real_part_product(coeffs[idx_coefs]);
        }
    }

    /* we now have rho and can evaluate lots of things. 
       -------- this is back to the original PACE code --------- */

#ifdef DEBUG_FORCES_CALCULATIONS
    printf("rhos = ");
    for(DENSITY_TYPE p =0; p<ndensity; ++p) printf(" %.20f ",rhos(p));
    printf("\n");
#endif


    // energy cutoff
    rho_cut = basis_set->rho_core_cutoffs(mu_i);
    drho_cut = basis_set->drho_core_cutoffs(mu_i);

    basis_set->inner_cutoff(rho_core, rho_cut, drho_cut, fcut, dfcut);
    basis_set->FS_values_and_derivatives(rhos, evdwl, dF_drho, ndensity);

    dF_drho_core = evdwl * dfcut + 1;
    for (DENSITY_TYPE p = 0; p < ndensity; ++p)
        dF_drho(p) *= fcut;
    evdwl_cut = evdwl * fcut + rho_core;

    // E0 shift 
    evdwl_cut += basis_set->E0vals(mu_i);

    /* I've moved this from below the weight calculation
          since I believe it only times the energy? the weights 
          are only needed for the forces?
          But I believe we could add a third timer for computing just 
          the weights; this will allow us to check better where the 
          bottleneck is.
    */
    energy_calc_timer.stop();

    forces_calc_loop_timer.start();


#ifdef DEBUG_FORCES_CALCULATIONS
    printf("dFrhos = ");
    for(DENSITY_TYPE p =0; p<ndensity; ++p) printf(" %f ",dF_drho(p));
    printf("\n");
#endif

    //ALGORITHM 3: Weights and theta calculation
    // rank = 1
    for (int f_ind = 0; f_ind < total_basis_size_rank1; ++f_ind) {
        ACECTildeBasisFunction *func = &basis_rank1[f_ind];
//        ndensity = func->ndensity;
        for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
            //for rank=1 (r=0) only 1 ms-combination exists (ms_ind=0), so index of func.ctildes is 0..ndensity-1
            weights_rank1(func->mus[0], func->ns[0] - 1) += dF_drho(p) * func->ctildes[p];
        }
    }

    /* --------- we now continue with the recursive code --------- */

    if (recursive) {
        /* STAGE 2:  BACKWARD PASS */
        int i1, i2;
        ACEComplex AA1{0.0};
        ACEComplex AA2{0.0};
        ACEComplex wcur{0.0};
        int num2_int = dag.get_num2_int();
        int num2_leaf = dag.get_num2_leaf();
        /* to prepare for the backward we first need to zero the weights */
        dag.w.fill({0.0});

        int * dag_nodes = dag.nodes.get_data();
        int idx_nodes = 2 * (num2_int + num2_leaf) - 1;

        DOUBLE_TYPE * dag_coefs = dag.coeffs.get_data();
        int idx_coefs = ndensity * (num2_int + num2_leaf) - 1;
        
        for (int idx = num1+num2_int+num2_leaf - 1; idx >= num1; idx--) {
            i2 = dag_nodes[idx_nodes]; idx_nodes--;
            i1 = dag_nodes[idx_nodes]; idx_nodes--;
            AA1 = dag.AAbuf(i1);
            AA2 = dag.AAbuf(i2);
            wcur = dag.w(idx);   // [***]
            for (int p = ndensity-1; p >= 0; p--, idx_coefs--) 
                wcur += dF_drho(p) * dag_coefs[idx_coefs];
            dag.w(i1) += wcur * AA2;   // TODO: replace with explicit muladd? 
            dag.w(i2) += wcur * AA1;
        }

        /*  [***]
         * Note that these weights don't really need to be stored for the 
         * leaf nodes. We tested splitting this for loop into two where 
         * for the leaf nodes the weight would just be initialized to 0.0
         * instead of reading from an array. The improvement was barely 
         * measurable, ca 3%, so we reverted to this simpler algorithm
         */


    } else {

        // non-recursive ACE.jl style implemenation of gradients, but with 
        // a backward differentiation approach to the prod-A 
        // (cf. Algorithm 3 in the manuscript)

        dag.w.fill({0.0});
        ACEComplex AAf{1.0}, AAb{1.0}, theta{0.0};        

        int *AAspec = jl_AAspec_flat.get_data();
        DOUBLE_TYPE *coeffs = jl_coeffs.get_data();
        int idx_spec = 0;
        int idx_coefs = 0;
        int order = 0;
        int max_order = jl_AAspec.get_dim(1);
        for (int iAA = 0; iAA < jl_AAspec.get_dim(0); iAA ++ ) {
            order = jl_orders(iAA);
            theta = 0.0; 
            for (int p = 0; p < ndensity; p++, idx_coefs++) 
                theta += dF_drho(p) * coeffs[idx_coefs];
            dA[0] = 1.0; 
            AAf = 1.0; 
            for (int t = 0; t < order-1; t++, idx_spec++) {
                spec[t] = AAspec[idx_spec]; 
                A_cache[t] = dag.AAbuf(spec[t]);
                AAf *= A_cache[t]; 
                dA[t+1] = AAf; 
            }
            spec[order-1] = AAspec[idx_spec]; idx_spec++;
            A_cache[order-1] = dag.AAbuf(spec[order-1]);
            AAb = 1.0;
            for (int t = order-1; t >= 1; t--) {
                AAb *= A_cache[t]; 
                dA[t-1] *= AAb; 
                dag.w(spec[t]) += theta * dA[t];
            }
            dag.w(spec[0]) += theta * dA[0];
        }

    }

    /* STAGE 3: 
     * get the gradients from the 1-particle basis gradients and write them 
     * into the dF/drho derivatives.
     */
    /* In order to reuse the original PACE code, we copy the weights back 
     * into the the PACE datastructure. */

    for (int idx = 0; idx < num1; idx++) {
        int m = dag.Aspec(idx, 3); 
        if (m >= 0) {
            weights(dag.Aspec(idx, 0),      // mu
                    dag.Aspec(idx, 1) - 1,  // n
                    dag.Aspec(idx, 2),      // l
                    m ) += dag.w(idx);
        } else {
            int factor = (m % 2 == 0 ? 1 : -1);
            weights(dag.Aspec(idx, 0),      // mu
                    dag.Aspec(idx, 1) - 1,  // n
                    dag.Aspec(idx, 2),      // l
                    -m ) += factor * dag.w(idx).conjugated();
        }
    }


    /* ------ From here we are now back to the original PACE code ---- */

// ==================== FORCES ====================
#ifdef PRINT_MAIN_STEPS
    printf("\nFORCE CALCULATION\n");
    printf("loop over neighbours\n");
#endif

// loop over neighbour atoms for force calculations
    for (jj = 0; jj < jnum_actual; ++jj) {
        mu_j = elements[jj];
        r_hat = rhats[jj];
        inv_r_norm = inv_r_norms[jj];

        Array2DLM<ACEComplex> &Y_cache_jj = Y_cache(jj);
        Array2DLM<ACEDYcomponent> &DY_cache_jj = DY_cache(jj);

#ifdef PRINT_LOOPS_INDICES
        printf("\nneighbour atom #%d\n", jj);
        printf("rhat = (%f, %f, %f)\n", r_hat[0], r_hat[1], r_hat[2]);
#endif

        forces_calc_neighbour_timer.start();

        f_ji[0] = f_ji[1] = f_ji[2] = 0;

//for rank = 1
        for (n = 0; n < nradbasei; ++n) {
            if (weights_rank1(mu_j, n) == 0)
                continue;
            auto &DG = DG_cache(jj, n);
            DGR = DG * Y00;
            DGR *= weights_rank1(mu_j, n);
#ifdef DEBUG_FORCES_CALCULATIONS
            printf("r=1: (n,l,m)=(%d, 0, 0)\n",n+1);
            printf("\tGR(n=%d, r=%f)=%f\n",n+1,r_norm, gr(n));
            printf("\tDGR(n=%d, r=%f)=%f\n",n+1,r_norm, dgr(n));
            printf("\tdF+=(%f, %f, %f)\n",DGR * r_hat[0], DGR * r_hat[1], DGR * r_hat[2]);
#endif
            f_ji[0] += DGR * r_hat[0];
            f_ji[1] += DGR * r_hat[1];
            f_ji[2] += DGR * r_hat[2];
        }

//for rank > 1
        for (n = 0; n < nradiali; n++) {
            for (l = 0; l <= lmaxi; l++) {
                R_over_r = R_cache(jj, n, l) * inv_r_norm;
                DR = DR_cache(jj, n, l);

                // for m>=0
                for (m = 0; m <= l; m++) { 
                    ACEComplex w = weights(mu_j, n, l, m);
                    if (w == 0)
                        continue;
                    //counting for -m cases if m>0
                    // if (m > 0) w *= 2;  // not needed for recursive eval

                    DY = DY_cache_jj(l, m);
                    Y_DR = Y_cache_jj(l, m) * DR;

                    grad_phi_nlm.a[0] = Y_DR * r_hat[0] + DY.a[0] * R_over_r;
                    grad_phi_nlm.a[1] = Y_DR * r_hat[1] + DY.a[1] * R_over_r;
                    grad_phi_nlm.a[2] = Y_DR * r_hat[2] + DY.a[2] * R_over_r;
#ifdef DEBUG_FORCES_CALCULATIONS
                    printf("d_phi(n=%d, l=%d, m=%d) = ((%f,%f), (%f,%f), (%f,%f))\n",n+1,l,m,
                           grad_phi_nlm.a[0].real, grad_phi_nlm.a[0].img,
                           grad_phi_nlm.a[1].real, grad_phi_nlm.a[1].img,
                           grad_phi_nlm.a[2].real, grad_phi_nlm.a[2].img);

                    printf("weights(n,l,m)(%d,%d,%d) = (%f,%f)\n", n+1, l, m, w.real, w.img);
                    //if (m>0) w*=2;
                    printf("dF(n,l,m)(%d, %d, %d) += (%f, %f, %f)\n", n + 1, l, m,
                           w.real_part_product(grad_phi_nlm.a[0]),
                           w.real_part_product(grad_phi_nlm.a[1]),
                           w.real_part_product(grad_phi_nlm.a[2])
                    );
#endif
// real-part multiplication only
                    f_ji[0] += w.real_part_product(grad_phi_nlm.a[0]);
                    f_ji[1] += w.real_part_product(grad_phi_nlm.a[1]);
                    f_ji[2] += w.real_part_product(grad_phi_nlm.a[2]);
                }
            }
        }


#ifdef PRINT_INTERMEDIATE_VALUES
        printf("f_ji(jj=%d, i=%d)=(%f, %f, %f)\n", jj, i,
               f_ji[0], f_ji[1], f_ji[2]
        );
#endif

        //hard-core repulsion
        DCR = DCR_cache(jj);
#ifdef   DEBUG_FORCES_CALCULATIONS
        printf("DCR = %f\n",DCR);
#endif
        f_ji[0] += dF_drho_core * DCR * r_hat[0];
        f_ji[1] += dF_drho_core * DCR * r_hat[1];
        f_ji[2] += dF_drho_core * DCR * r_hat[2];
#ifdef PRINT_INTERMEDIATE_VALUES
        printf("with core-repulsion\n");
        printf("f_ji(jj=%d, i=%d)=(%f, %f, %f)\n", jj, i,
               f_ji[0], f_ji[1], f_ji[2]
        );
        printf("neighbour_index_mapping[jj=%d]=%d\n",jj,neighbour_index_mapping[jj]);
#endif

        neighbours_forces(neighbour_index_mapping[jj], 0) = f_ji[0];
        neighbours_forces(neighbour_index_mapping[jj], 1) = f_ji[1];
        neighbours_forces(neighbour_index_mapping[jj], 2) = f_ji[2];

        forces_calc_neighbour_timer.stop();
    }// end loop over neighbour atoms for forces

    forces_calc_loop_timer.stop();

    //now, energies and forces are ready
    //energies(i) = evdwl + rho_core;
    e_atom = evdwl_cut;

#ifdef PRINT_INTERMEDIATE_VALUES
    printf("energies(i) = FS(...rho_p_accum...) = %f\n", evdwl);
#endif
    per_atom_calc_timer.stop();
}