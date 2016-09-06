#pragma once
#ifndef __CHOL_U_HPP__
#define __CHOL_U_HPP__

/// \file chol_u.hpp
/// \brief Upper Cholesky factorization variations
/// \author Kyungjoo Kim (kyukim@sandia.gov)

// testing task-data parallelism
// #include "chol_u_unblocked_dummy.hpp"

// flame style implementation
//#include "chol_unblocked.hpp"  
//#include "chol_u_blocked.hpp"

// triple for loop
#include "chol_u_unblocked_opt1.hpp"
#include "chol_u_unblocked_opt2.hpp"

// partitioned block algorithms: see control.hpp
#include "chol_u_right_look_by_blocks.hpp"

#endif
