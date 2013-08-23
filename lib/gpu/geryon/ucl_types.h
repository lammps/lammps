/***************************************************************************
                                 ucl_types.h
                             -------------------
                               W. Michael Brown

  Data type definitions for Coprocessor library

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Mon Jan 4 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifndef UCL_TYPES_H
#define UCL_TYPES_H

// Assign an integer id based on the data type: (int, float, double, etc)
template <class eltype> struct _UCL_DATA_ID;
template <> struct _UCL_DATA_ID<double> { 
  enum { id=1 };
  static inline const char * name() { return "double"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=double"; }  
};
template <> struct _UCL_DATA_ID<float> { 
  enum { id=2 };
  static inline const char * name() { return "float"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=float"; }  
};
template <> struct _UCL_DATA_ID<unsigned> { 
  enum { id=3 };
  static inline const char * name() { return "unsigned"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=unsigned"; }  
};
template <> struct _UCL_DATA_ID<int> { 
  enum { id=4 };
  static inline const char * name() { return "int"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=int"; }  
};
template <> struct _UCL_DATA_ID<char> { 
  enum { id=5 };
  static inline const char * name() { return "char"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=char"; }  
};
template <> struct _UCL_DATA_ID<unsigned char> { 
  enum { id=6 };
  static inline const char * name() { return "unsigned char"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=unsigned char"; }  
};
template <> struct _UCL_DATA_ID<short> { 
  enum { id=7 };
  static inline const char * name() { return "short"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=short"; }  
};
template <> struct _UCL_DATA_ID<unsigned short> { 
  enum { id=8 };
  static inline const char * name() { return "unsigned short"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=unsigned short"; }  
};
template <> struct _UCL_DATA_ID<long> { 
  enum { id=9 };
  static inline const char * name() { return "long"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=long"; }  
};
template <> struct _UCL_DATA_ID<unsigned long> { 
  enum { id=10 };
  static inline const char * name() { return "unsigned long"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=unsigned long"; }  
};
template <> struct _UCL_DATA_ID<long double> { 
  enum { id=11 };
  static inline const char * name() { return "long double"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=long double"; }  
};
template <class eltype> struct _UCL_DATA_ID { 
  enum { id=0 };
  static inline const char * name() { return "error_type"; }  
  static inline const char * numtyp_flag() { return "-D NUMTYP=error_type"; }  
};

// Host memory allocation types
enum UCL_MEMOPT {
  UCL_WRITE_ONLY,     ///< Allow any optimizations for memory that is write only
  UCL_READ_ONLY,      ///< Allow any optimizations for memory that is read only
  UCL_READ_WRITE,     ///< Allow read and write
  UCL_NOT_PINNED,     ///< Host memory is not to be pinned
  UCL_VIEW,           ///< View of another memory allocation
  UCL_NOT_SPECIFIED
};

enum UCL_DEVICE_TYPE { 
  UCL_DEFAULT,        ///< Unknown device type
  UCL_CPU,            ///< Device is a CPU
  UCL_GPU,            ///< Device is a GPU
  UCL_ACCELERATOR     ///< Device is an Accelerator
};

enum UCL_ERROR_FLAG {
  UCL_SUCCESS,            ///< No error
  UCL_ERROR,              ///< Unqualified error
  UCL_FILE_NOT_FOUND,     ///< File not found
  UCL_FUNCTION_NOT_FOUND, ///< Kernel function not found
  UCL_COMPILE_ERROR,      ///< Error compiling kernel
  UCL_MEMORY_ERROR
};  

template <class numtyp>
const char * ucl_template_name() { return _UCL_DATA_ID<numtyp>::name(); }

template <class t1, class t2> struct ucl_same_type;

template <> struct ucl_same_type<bool,bool> { enum { ans=1 }; };
template <> struct ucl_same_type<char,char> { enum { ans=1 }; };
template <> struct ucl_same_type<unsigned char,unsigned char> { enum { ans=1 }; };
template <> struct ucl_same_type<int,int> { enum { ans=1 }; };
template <> struct ucl_same_type<unsigned,unsigned> { enum { ans=1 }; };
template <> struct ucl_same_type<short,short> { enum { ans=1 }; };
template <> struct ucl_same_type<unsigned short,unsigned short> { enum { ans=1 }; };
template <> struct ucl_same_type<long,long> { enum { ans=1 }; };
template <> struct ucl_same_type<unsigned long,unsigned long> { enum { ans=1 }; };
template <> struct ucl_same_type<float,float> { enum { ans=1 }; };
template <> struct ucl_same_type<double,double> { enum { ans=1 }; };
template <> struct ucl_same_type<long double,long double> { enum { ans=1 }; };

template <> struct ucl_same_type<const bool,bool> { enum { ans=1 }; };
template <> struct ucl_same_type<const char,char> { enum { ans=1 }; };
template <> struct ucl_same_type<const unsigned char,unsigned char> { enum { ans=1 }; };
template <> struct ucl_same_type<const int,int> { enum { ans=1 }; };
template <> struct ucl_same_type<const unsigned,unsigned> { enum { ans=1 }; };
template <> struct ucl_same_type<const short,short> { enum { ans=1 }; };
template <> struct ucl_same_type<const unsigned short,unsigned short> { enum { ans=1 }; };
template <> struct ucl_same_type<const long,long> { enum { ans=1 }; };
template <> struct ucl_same_type<const unsigned long,unsigned long> { enum { ans=1 }; };
template <> struct ucl_same_type<const float,float> { enum { ans=1 }; };
template <> struct ucl_same_type<const double,double> { enum { ans=1 }; };
template <> struct ucl_same_type<const long double,long double> { enum { ans=1 }; };

template <> struct ucl_same_type<bool,const bool> { enum { ans=1 }; };
template <> struct ucl_same_type<char,const char> { enum { ans=1 }; };
template <> struct ucl_same_type<unsigned char,const unsigned char> { enum { ans=1 }; };
template <> struct ucl_same_type<int,const int> { enum { ans=1 }; };
template <> struct ucl_same_type<unsigned,const unsigned> { enum { ans=1 }; };
template <> struct ucl_same_type<short,const short> { enum { ans=1 }; };
template <> struct ucl_same_type<unsigned short,const unsigned short> { enum { ans=1 }; };
template <> struct ucl_same_type<long,const long> { enum { ans=1 }; };
template <> struct ucl_same_type<unsigned long,const unsigned long> { enum { ans=1 }; };
template <> struct ucl_same_type<float,const float> { enum { ans=1 }; };
template <> struct ucl_same_type<double,const double> { enum { ans=1 }; };
template <> struct ucl_same_type<long double,const long double> { enum { ans=1 }; };

template <> struct ucl_same_type<const bool,const bool> { enum { ans=1 }; };
template <> struct ucl_same_type<const char,const char> { enum { ans=1 }; };
template <> struct ucl_same_type<const unsigned char,const unsigned char> { enum { ans=1 }; };
template <> struct ucl_same_type<const int,const int> { enum { ans=1 }; };
template <> struct ucl_same_type<const unsigned,const unsigned> { enum { ans=1 }; };
template <> struct ucl_same_type<const short,const short> { enum { ans=1 }; };
template <> struct ucl_same_type<const unsigned short,const unsigned short> { enum { ans=1 }; };
template <> struct ucl_same_type<const long,const long> { enum { ans=1 }; };
template <> struct ucl_same_type<const unsigned long,const unsigned long> { enum { ans=1 }; };
template <> struct ucl_same_type<const float,const float> { enum { ans=1 }; };
template <> struct ucl_same_type<const double,const double> { enum { ans=1 }; };
template <> struct ucl_same_type<const long double,const long double> { enum { ans=1 }; };

template <class t1, class t2> struct ucl_same_type { enum { ans=0 }; };

#endif

