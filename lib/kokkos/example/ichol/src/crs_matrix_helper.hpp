#pragma once
#ifndef __CRS_MATRIX_HELPER_HPP__
#define __CRS_MATRIX_HELPER_HPP__

/// \file crs_matrix_helper.hpp
/// \brief This file includes utility functions to convert between flat and hierarchical matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)  

#include "util.hpp"

namespace Tacho { 

  using namespace std;

  class CrsMatrixHelper {
  public:

    template< typename CrsHierBase >
    static int fillRowViewArray( CrsHierBase & HU );

    template<typename CrsFlatBase>
    static int
    filterZeros(CrsFlatBase &flat);
    
    /// \brief Transform a scalar flat matrix to hierarchical matrix of matrices 1x1; testing only.
    template<typename CrsFlatBase,
             typename CrsHierBase>
    static int
    flat2hier(CrsFlatBase &flat, 
              CrsHierBase &hier);

    /// \brief Transform a scalar flat matrix to upper hierarchical matrix given scotch info. 
    template<typename CrsFlatBase,
             typename CrsHierBase,
             typename HostOrdinalTypeArray >
    static int
    flat2hier(int uplo, 
              CrsFlatBase &flat, 
              CrsHierBase &hier,
              const typename CrsHierBase::ordinal_type       nblks,
              const HostOrdinalTypeArray range,
              const HostOrdinalTypeArray tree);

    /// \brief Transform a scalar flat matrix to upper hierarchical matrix given scotch info. 
    template<typename CrsFlatBase,
             typename CrsHierBase,
             typename HostOrdinalTypeArray >
    static int
    flat2hier_upper(CrsFlatBase &flat, 
                    CrsHierBase &hier,
                    const typename CrsHierBase::ordinal_type       nblks,
                    const HostOrdinalTypeArray range,
                    const HostOrdinalTypeArray tree);

    /// \brief Transform a scalar flat matrix to lower hierarchical matrix given scotch info. 
    template<typename CrsFlatBase,
             typename CrsHierBase,
             typename HostOrdinalTypeArray >
    static int
    flat2hier_lower(CrsFlatBase &flat, 
                    CrsHierBase &hier,
                    const typename CrsHierBase::ordinal_type       nblks,
                    const HostOrdinalTypeArray range,
                    const HostOrdinalTypeArray tree);
  };

}

#include "crs_matrix_helper_impl.hpp"

#endif
