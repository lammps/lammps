#pragma once
#ifndef __COO_HPP__
#define __COO_HPP__

/// \file coo.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 
  
  using namespace std;

  /// \class Coo
  /// \brief Sparse coordinate format; (i, j, val).
  template<typename CrsMatType>
  class Coo {
  public:
    typedef typename CrsMatType::ordinal_type ordinal_type;
    typedef typename CrsMatType::value_type   value_type;

  public:
    ordinal_type _i,_j;
    value_type _val;

  public:
    ordinal_type& Row() { return _i;   } 
    ordinal_type& Col() { return _j;   }
    value_type&   Val() { return _val; }

    ordinal_type  Row() const { return _i;   } 
    ordinal_type  Col() const { return _j;   }
    value_type    Val() const { return _val; }
    
    Coo() {}

    Coo(const ordinal_type i, 
        const ordinal_type j, 
        const value_type val) 
      : _i(i),
        _j(j),
        _val(val) 
    { }

    Coo(const Coo& b)
      : _i(b._i),
        _j(b._j),
        _val(b._val) 
    { }

    Coo<CrsMatType>& operator=(const Coo<CrsMatType> &y) {
      this->_i = y._i;
      this->_j = y._j;
      this->_val = y._val;

      return *this;
    }

    /// \brief Compare "less" index i and j only.
    bool operator<(const Coo<CrsMatType> &y) const {
      ordinal_type r_val = (this->_i - y._i);
      return (r_val == 0 ? this->_j < y._j : r_val < 0);
    }  
    
    /// \brief Compare "equality" only index i and j.
    bool operator==(const Coo<CrsMatType> &y) const {
      return (this->_i == y._i) && (this->_j == y._j);
    }  
 
    /// \brief Compare "in-equality" only index i and j.   
    bool operator!=(const Coo<CrsMatType> &y) const {
      return !(*this == y);
    }  
  };
  
}
#endif
