#pragma once
#ifndef __CRS_MATRIX_BASE_IMPL_HPP__
#define __CRS_MATRIX_BASE_IMPL_HPP__

/// \file crs_matrix_base_impl.hpp
/// \brief Implementation of external interfaces to CrsMatrixBase
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 

  using namespace std;

  template<typename VT,
           typename OT,
           typename ST,
           typename SpT,
           typename MT>
  inline int 
  CrsMatrixBase<VT,OT,ST,SpT,MT>::importMatrixMarket(ifstream &file) {
    // skip initial title comments
    {
      ordinal_type m, n;
      size_type nnz;
          
      while (file.good()) {
        char c = file.peek();
        if (c == '%' || c == '\n') {
          file.ignore(256, '\n');
          continue;
        }
        break;
      }
          
      // read matrix specification
      file >> m >> n >> nnz;
          
      // construct workspace and set variables
      createInternalArrays(m, n, nnz);
    }

    // read the coordinate format (matrix-market)
    vector<ijv_type> mm; 
    mm.reserve(_nnz);
    {
      // matrix market use one base index
      const ordinal_type mm_base = 1; 

      for (size_type i=0;i<_nnz;++i) {
        ijv_type aij;
        file >> aij.Row() >> aij.Col() >> aij.Val();

        // one base to zero base
        aij.Row() -= mm_base;
        aij.Col() -= mm_base;
            
        mm.push_back(aij);
      }
      sort(mm.begin(), mm.end(), less<ijv_type>());
    }

    // change mm to crs
    {
      ordinal_type ii = 0;
      size_type jj = 0;

      ijv_type prev = mm[0];
      _ap[ii++] = 0;
      _aj[jj] = prev.Col();
      _ax[jj] = prev.Val();
      ++jj;

      for (typename vector<ijv_type>::iterator it=(mm.begin()+1);it<mm.end();++it) {
        ijv_type aij = (*it);
        
        // row index
        if (aij.Row() != prev.Row()) {
          _ap[ii++] = jj; 
        }
            
        if (aij == prev) {
          --jj;
          _aj[jj]  = aij.Col();
          _ax[jj] += aij.Val();
        } else {
          _aj[jj] = aij.Col();
          _ax[jj] = aij.Val();
        }
        ++jj;
          
        prev = aij;
      }
          
      // add the last index to terminate the storage
      _ap[ii++] = jj;
      _nnz = jj;
    }
      
    return 0;
  }
  
}


#endif
