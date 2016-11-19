#pragma once
#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <stdio.h>
#include <string.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <memory>

#include <cmath>
#include <complex>

#include <limits>

/// \file util.hpp
/// \brief Utility functions and constant integer class like an enum class.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
/// This provides utility functions for implementing mini-app for incomplete
/// sparse matrix factorization with task-data parallelism e.g., parameter
/// classes, error handling, ostream << overloading.
///
/// Note: The reference of the "static const int" members in the enum-like
/// classes should not be used as function arguments but their values only.


using namespace std;

namespace Tacho {

#undef CHKERR
#define CHKERR(ierr)                                                    \
  if (ierr != 0) { cout << endl << ">> Error in " << __FILE__ << ", " << __LINE__ << " : " << ierr << endl; }

#define MSG_NOT_YET_IMPLEMENTED ">> Not yet implemented"
#define MSG_INVALID_INPUT(what) ">> Invaid input argument: " #what
#define MSG_INVALID_TEMPLATE_ARGS ">> Invaid template arguments"
#define ERROR(msg)                                                      \
  { cout << endl << ">> Error in " << __FILE__ << ", " << __LINE__ << endl << msg << endl; }

  // control id
#undef  Ctrl
#define Ctrl(name,algo,variant) name<algo,variant>

  // control leaf
#undef CtrlComponent
#define CtrlComponent(name,algo,variant,component,id)                  \
  Ctrl(name,algo,variant)::component[id]

  // control recursion
#undef CtrlDetail
#define CtrlDetail(name,algo,variant,component) \
  CtrlComponent(name,algo,variant,component,0),CtrlComponent(name,algo,variant,component,1),name

  /// \class GraphHelper
  class GraphHelper {
  public:
    static const int DefaultRandomSeed = -1;
  };


  /// \class Partition
  /// \brief Matrix partition parameters.
  class Partition {
  public:
    static const int Top         = 101;
    static const int Bottom      = 102;

    static const int Left        = 201;
    static const int Right       = 202;

    static const int TopLeft     = 401;
    static const int TopRight    = 402;
    static const int BottomLeft  = 403;
    static const int BottomRight = 404;
  };

  /// \class Uplo
  /// \brief Matrix upper/lower parameters.
  class Uplo {
  public:
    static const int Upper = 501;
    static const int Lower = 502;
  };

  /// \class Side
  /// \brief Matrix left/right parameters.
  class Side {
  public:
    static const int Left  = 601;
    static const int Right = 602;
  };

  /// \class Diag
  /// \brief Matrix unit/non-unit diag parameters.
  class Diag {
  public:
    static const int Unit    = 701;
    static const int NonUnit = 702;
  };

  /// \class Trans
  /// \brief Matrix upper/lower parameters.
  class Trans {
  public:
    static const int Transpose     = 801;
    static const int ConjTranspose = 802;
    static const int NoTranspose   = 803;
  };

  /// \class Loop
  /// \brief outer/innner parameters
  class Loop {
  public:
    static const int Outer = 901;
    static const int Inner = 902;
    static const int Fused = 903;
  };

  class Variant {
  public:
    static const int One   = 1;
    static const int Two   = 2;
    static const int Three = 3;
    static const int Four  = 4;
  };

  /// \class AlgoChol
  /// \brief Algorithmic variants in sparse factorization and sparse BLAS operations.
  class AlgoChol {
  public:
    // One side factorization on flat matrices
    static const int Dummy                  = 1000;
    static const int Unblocked              = 1001;
    static const int UnblockedOpt           = 1002;
    static const int Blocked                = 1101; // testing only

    static const int RightLookByBlocks      = 1201; // backbone structure is right looking
    static const int ByBlocks               = RightLookByBlocks;

    static const int NestedDenseBlock       = 1211;
    static const int NestedDenseByBlocks    = 1212;

    static const int RightLookDenseByBlocks = 1221;
    static const int DenseByBlocks          = RightLookDenseByBlocks;

    static const int ExternalLapack         = 1231;
    static const int ExternalPardiso        = 1232;
  };

  // aliasing name space
  typedef AlgoChol AlgoTriSolve;

  class AlgoBlasLeaf {
  public:
    // One side factorization on flat matrices
    static const int ForFactorBlocked = 2001;

    // B and C are dense matrices and used for solve phase
    static const int ForTriSolveBlocked = 2011;

    static const int ExternalBlas = 2021;
  };

  class AlgoGemm : public AlgoBlasLeaf {
  public:
    static const int DenseByBlocks = 2101;
  };

  class AlgoTrsm : public AlgoBlasLeaf {
  public:
    static const int DenseByBlocks = 2201;
  };

  class AlgoHerk : public AlgoBlasLeaf {
  public:
    static const int DenseByBlocks = 2301;
  };

  /// \brief Interface for overloaded stream operators.
  template<typename T>
  inline
  ostream& operator<<(ostream &os, const unique_ptr<T> &p) {
    return p->showMe(os);
  }

  /// \class Disp
  /// \brief Interface for the stream operator.
  class Disp {
    friend ostream& operator<<(ostream &os, const Disp &disp);
  public:
    Disp() { }
    virtual ostream& showMe(ostream &os) const {
      return os;
    }
  };

  /// \brief Implementation of the overloaded stream operator.
  inline
  ostream& operator<<(ostream &os, const Disp &disp) {
    return disp.showMe(os);
  }

  template<typename T> struct NumericTraits {};

  template<>
  struct NumericTraits<float> {
    typedef float real_type;
    static real_type epsilon() { return numeric_limits<float>::epsilon(); }
  };
  template<>
  struct NumericTraits<double> {
    typedef double real_type;
    static real_type epsilon() { return numeric_limits<double>::epsilon(); }
  };
  template<>
  struct NumericTraits<complex<float> > {
    typedef float real_type;
    static real_type epsilon() { return numeric_limits<float>::epsilon(); }
  };
  template<>
  struct NumericTraits<complex<double> > {
    typedef double real_type;
    static real_type epsilon() { return numeric_limits<double>::epsilon(); }
  };

}

#endif
