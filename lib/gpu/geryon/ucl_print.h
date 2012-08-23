/***************************************************************************
                                 ucl_print.h
                             -------------------
                               W. Michael Brown

  Routines for printing debugging output for matrix/vector data

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Mon Jan 11 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the Simplified BSD License.
   ----------------------------------------------------------------------- */
   
// Only allow this file to be included by nvc_memory.h and ocl_memory.h
#ifdef UCL_PRINT_ALLOW

template <int mem> struct _ucl_print;
template <> struct _ucl_print<1> {
  template <class mat_type>
  static inline void p(mat_type &mat, const size_t n, std::ostream &out,
                       const std::string delim) {
    for (size_t i=0; i<n-1; i++)
      out << mat[i] << delim;
    out << mat[n-1];
  }
  template <class mat_type>
  static inline void p(const mat_type &mat, const size_t n, std::ostream &out,
                       const std::string delim, UCL_Device &dev) {
    p(mat,n,out,delim);
  }
  template <class mat_type>
  static inline void p(mat_type &mat, const size_t rows, const size_t cols,
                       std::ostream &out, const std::string delim, 
                       const std::string row_delim) {
    int offset=0;
    int row_size=cols;
    if (mat_type::VECTOR==0)
      row_size=mat.row_size();
    for (size_t j=0; j<rows; j++) {
      size_t lend=offset+cols-1;
      for (size_t i=offset; i<lend; i++)
        out << mat[i] << delim;
      out << mat[lend];
      if (j!=rows-1)
        out << row_delim;
      offset+=row_size;
    }
  }
  template <class mat_type>
  static inline void p(const mat_type &mat,const size_t rows,const size_t cols,
                       std::ostream &out,const std::string delim, 
                       const std::string row_delim, UCL_Device &dev) {
    p(mat,rows,cols,out,delim,row_delim);                       
  }
};
      
template <int mem> struct _ucl_print {
  template <class mat_type>
  static inline void p(mat_type &mat, const size_t n, std::ostream &out,
                       const std::string delim) {
    UCL_H_Vec<typename mat_type::data_type> temp;
    temp.alloc(n,mat);
    ucl_copy(temp,mat,n,false);
    _ucl_print<1>::p(temp,n,out,delim);
  }
  template <class mat_type>
  static inline void p(const mat_type &mat, const size_t n, std::ostream &out,
                       const std::string delim, UCL_Device &dev) {
    UCL_H_Vec<typename mat_type::data_type> temp;
    temp.alloc(n,dev);
    ucl_copy(temp,mat,n,false);
    _ucl_print<1>::p(temp,n,out,delim);
  }
  template <class mat_type>
  static inline void p(mat_type &mat, const size_t rows, const size_t cols,
                       std::ostream &out, const std::string delim, 
                       const std::string row_delim) {
    UCL_H_Vec<typename mat_type::data_type> temp;
    temp.alloc(mat.rows()*mat.cols(),mat);
    if (mat_type::VECTOR==1)
      ucl_copy(temp,mat,rows*cols,false);
    else
      ucl_copy(temp,mat,rows,cols,false);
    _ucl_print<1>::p(temp,rows,cols,out,delim,row_delim);      
  }
  template <class mat_type>
  static inline void p(const mat_type &mat, const size_t rows, 
                       const size_t cols,std::ostream &out,
                       const std::string delim, 
                       const std::string row_delim, UCL_Device &dev) {
    UCL_H_Vec<typename mat_type::data_type> temp;
    temp.alloc(mat.rows()*mat.cols(),dev);
    if (mat_type::VECTOR==1)
      ucl_copy(temp,mat,rows*cols,false);
    else
      ucl_copy(temp,mat,rows,cols,false);
    _ucl_print<1>::p(temp,rows,cols,out,delim,row_delim);      
  }
};                   

// -------------------------------------------------------------------------
// - Non-const routines that do not require a device object
// -------------------------------------------------------------------------

/// Outputs n elements of mat delimited by the string delim
template <class mat_type>
inline void ucl_print(mat_type &mat, const size_t n, std::ostream &out,
                      const std::string delim) {
  if (n>mat.numel()) {
    std::cerr << "Attempted to ucl_print " << n << " elements of matrix "
              << "that only has " << mat.numel() << " elements.";
    UCL_GERYON_EXIT;
  }
  _ucl_print<mat_type::MEM_TYPE>::p(mat,n,out,delim);
}
  
/// Outputs n elements of mat delimited by a space
template <class mat_type>
inline void ucl_print(mat_type &mat, const size_t n, std::ostream &out) {
  ucl_print(mat,n,out," ");
}
  
/// Outputs n elements of mat delimited by a space to standard out
template <class mat_type>
inline void ucl_print(mat_type &mat, const size_t n) {
  ucl_print(mat,n,std::cout," ");
}

/// Outputs upper left rows and cols of mat delimited by the string delim
template <class mat_type>
inline void ucl_print(mat_type &mat, const size_t rows, const size_t cols,
                      std::ostream &out, const std::string delim, 
                      const std::string row_delim) {                      
  if (rows*cols>mat.numel()) {
    std::cerr << "Attempted to ucl_print " << rows*cols << " elements of matrix "
              << "that only has " << mat.numel() << " elements.";
    UCL_GERYON_EXIT;
  }
  _ucl_print<mat_type::MEM_TYPE>::p(mat,rows,cols,out,delim,row_delim);
}
  
/// Outputs upper left rows and cols of mat delimited by a space
template <class mat_type>
inline void ucl_print(mat_type &mat, const size_t rows, const size_t cols,
                      std::ostream &out) {
  ucl_print(mat,rows,cols,out," ","\n");
}
  
/// Outputs  upper left rows and cols of mat delimited by a space to std out
template <class mat_type>
inline void ucl_print(mat_type &mat, const size_t rows, 
                      const size_t cols) {
  ucl_print(mat,rows,cols,std::cout," ","\n");
}

/// Outputs mat delimited by a space to standard out
template <class mat_type>
inline void ucl_print(mat_type &mat) {
  ucl_print(mat,std::cout);
}

/// Outputs mat delimited by a space
template <class mat_type>
inline void ucl_print(mat_type &mat, std::ostream &out) {
  if (mat_type::VECTOR==1)
    ucl_print(mat,mat.cols(),out," ");
  else
    ucl_print(mat,mat.rows(),mat.cols(),out," ","\n");
}
  
// -------------------------------------------------------------------------
// - Const routines that do not require a device object
// -------------------------------------------------------------------------

/// Outputs n elements of mat delimited by the string delim
template <class mat_type>
inline void ucl_print(const mat_type &mat, const size_t n, std::ostream &out,
                      const std::string delim, UCL_Device &dev) {
  if (n>mat.numel()) {
    std::cerr << "Attempted to ucl_print " << n << " elements of matrix "
              << "that only has " << mat.numel() << " elements.";
    UCL_GERYON_EXIT;
  }
  _ucl_print<mat_type::MEM_TYPE>::p(mat,n,out,delim,dev);
}
  
/// Outputs n elements of mat delimited by a space
template <class mat_type>
inline void ucl_print(const mat_type &mat, const size_t n, std::ostream &out,
                      UCL_Device &dev) {
  ucl_print(mat,n,out," ",dev);
}
  
/// Outputs n elements of mat delimited by a space to standard out
template <class mat_type>
inline void ucl_print(const mat_type &mat, const size_t n,
                      UCL_Device &dev) {
  ucl_print(mat,n,std::cout," ",dev);
}

/// Outputs upper left rows and cols of mat delimited by the string delim
template <class mat_type>
inline void ucl_print(const mat_type &mat,const size_t rows,const size_t cols,
                      std::ostream &out, const std::string delim, 
                      const std::string row_delim, UCL_Device &dev) {
  if (rows*cols>mat.numel()) {
    std::cerr << "Attempted to ucl_print " << rows*cols << " elements of matrix "
              << "that only has " << mat.numel() << " elements.";
    UCL_GERYON_EXIT;
  }
  _ucl_print<mat_type::MEM_TYPE>::p(mat,rows,cols,out,delim,row_delim,dev);
}
  
/// Outputs upper left rows and cols of mat delimited by a space
template <class mat_type>
inline void ucl_print(const mat_type &mat,const size_t rows,const size_t cols,
                      std::ostream &out, UCL_Device &dev) {
  ucl_print(mat,rows,cols,out," ","\n",dev);
}
  
/// Outputs  upper left rows and cols of mat delimited by a space to std out
template <class mat_type>
inline void ucl_print(const mat_type &mat, const size_t rows, 
                      const size_t cols, UCL_Device &dev) {
  ucl_print(mat,rows,cols,std::cout," ","\n",dev);
}

/// Outputs mat delimited by a space to standard out
template <class mat_type>
inline void ucl_print(const mat_type &mat, UCL_Device &dev) {
  ucl_print(mat,std::cout,dev);
}

/// Outputs mat delimited by a space
template <class mat_type>
inline void ucl_print(const mat_type &mat, std::ostream &out, UCL_Device &dev) {
  if (mat_type::VECTOR==1)
    ucl_print(mat,mat.cols(),out," ",dev);
  else
    ucl_print(mat,mat.rows(),mat.cols(),out," ","\n",dev);
}

// -------------------------------------------------------------------------
// - Operator << Overloading
// -------------------------------------------------------------------------

template <class numtyp>
inline std::ostream & operator << (std::ostream &out, UCL_H_Vec<numtyp> &mat)
  { ucl_print(mat,out); return out; } 

template <class numtyp>
inline std::ostream & operator << (std::ostream &out, UCL_H_Mat<numtyp> &mat)
  { ucl_print(mat,out); return out; } 

template <class numtyp>
inline std::ostream & operator << (std::ostream &out, UCL_D_Vec<numtyp> &mat)
  { ucl_print(mat,out); return out; } 

template <class numtyp>
inline std::ostream & operator << (std::ostream &out, UCL_D_Mat<numtyp> &mat)
  { ucl_print(mat,out); return out; } 


template <class t1, class t2>
inline std::ostream & operator << (std::ostream &out, UCL_Vector<t1,t2> &mat)
  { ucl_print(mat.host,out); return out; } 

template <class t1, class t2>
inline std::ostream & operator << (std::ostream &out, UCL_Matrix<t1,t2> &mat)
  { ucl_print(mat.host,out); return out; } 

#endif
