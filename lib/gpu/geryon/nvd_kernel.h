/***************************************************************************
                                nvd_kernel.h
                             -------------------
                               W. Michael Brown

  Utilities for dealing with CUDA Driver kernels

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Tue Feb 9 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifndef NVD_KERNEL
#define NVD_KERNEL

#include "nvd_device.h"
#include <fstream>
#include <cstdio>

namespace ucl_cudadr {

class UCL_Texture;
template <class numtyp> class UCL_D_Vec;
template <class numtyp> class UCL_D_Mat;
template <class hosttype, class devtype> class UCL_Vector;
template <class hosttype, class devtype> class UCL_Matrix;
#define UCL_MAX_KERNEL_ARGS 256

/// Class storing 1 or more kernel functions from a single string or file
class UCL_Program {
 public:
  inline UCL_Program(UCL_Device &device) { _cq=device.cq(); }
  inline UCL_Program(UCL_Device &device, const void *program,
                     const char *flags="", std::string *log=nullptr) {
    _cq=device.cq();
    init(device);
    load_string(program,flags,log);
  }

  inline ~UCL_Program() {}

  /// Initialize the program with a device
  inline void init(UCL_Device &device) { _cq=device.cq(); }

  /// Clear any data associated with program
  /** \note Must call init() after each clear **/
  inline void clear() { }

  /// Load a program from a file and compile with flags
  inline int load(const char *filename, const char *flags="",
                  std::string *log=nullptr) {
    std::ifstream in(filename);
    if (!in || in.is_open()==false) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not open kernel file: "
                << filename << std::endl;
      UCL_GERYON_EXIT;
      #endif
      return UCL_FILE_NOT_FOUND;
    }

    std::string program((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    in.close();
    return load_string(program.c_str(),flags,log);
  }

  /// Load a program from a string and compile with flags
  inline int load_string(const void *program, const char *flags="",
                         std::string *log=nullptr, FILE* foutput=nullptr) {
    if (std::string(flags)=="BINARY")
      return load_binary((const char *)program);
    const unsigned int num_opts=2;
    CUjit_option options[num_opts];
    void *values[num_opts];

    // set up size of compilation log buffer
    options[0] = CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES;
    values[0] = (void *)(int)10240;
    // set up pointer to the compilation log buffer
    options[1] = CU_JIT_INFO_LOG_BUFFER;
    char clog[10240];
    values[1] = clog;

    CUresult err=cuModuleLoadDataEx(&_module,program,num_opts,
                                    options,(void **)values);

    if (log!=nullptr)
      *log=std::string(clog);

    if (err != CUDA_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------------\n"
                << " UCL Error: Error compiling PTX Program...\n"
                << "----------------------------------------------------------\n";
      std::cerr << log << std::endl
                << "----------------------------------------------------------\n\n";
      #endif
      if (foutput != nullptr) {
        fprintf(foutput,"\n\n");
        fprintf(foutput, "----------------------------------------------------------\n");
        fprintf(foutput, " UCL Error: Error compiling PTX Program...\n");
        fprintf(foutput, "----------------------------------------------------------\n");
        fprintf(foutput, "%s\n",log->c_str());
        fprintf(foutput, "----------------------------------------------------------\n");
        fprintf(foutput,"\n\n");
      }
      return UCL_COMPILE_ERROR;
    }

    return UCL_SUCCESS;
  }

  /// Load a precompiled program from a file
  inline int load_binary(const char *filename) {
    CUmodule _module;
    CUresult err = cuModuleLoad(&_module,filename);
    if (err==301) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not open binary kernel file: "
                << filename << std::endl;
      UCL_GERYON_EXIT;
      #endif
      return UCL_FILE_NOT_FOUND;
    } else if (err!=CUDA_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Error loading binary kernel file: "
                << filename << std::endl;
      UCL_GERYON_EXIT;
      #endif
      return UCL_FILE_NOT_FOUND;
    }
    //int ucl_error=UCL_SUCCESS;
    //if (err==301)
    //  return UCL_FILE_NOT_FOUND;
    //else if (err!=CUDA_SUCCESS)
    //  return UCL_ERROR;
    return UCL_SUCCESS;
  }

  /// Return the default command queue/stream associated with this data
  inline command_queue & cq() { return _cq; }

  friend class UCL_Kernel;
 private:
  CUmodule _module;
  CUstream _cq;
  friend class UCL_Texture;
  friend class UCL_Const;
};

/// Class for dealing with CUDA Driver kernels
class UCL_Kernel {
 public:
  UCL_Kernel() : _dimensions(1), _num_args(0) {
    #if CUDA_VERSION < 4000
    _param_size=0;
    #endif
    _num_blocks[0]=0;
  }

  UCL_Kernel(UCL_Program &program, const char *function) :
    _dimensions(1), _num_args(0) {
    #if CUDA_VERSION < 4000
    _param_size=0;
    #endif
    _num_blocks[0]=0;
    set_function(program,function);
    _cq=program._cq;
  }

  ~UCL_Kernel() {}

  /// Clear any function associated with the kernel
  inline void clear() { }

  /// Get the kernel function from a program
  /** \ret UCL_ERROR_FLAG (UCL_SUCCESS, UCL_FILE_NOT_FOUND, UCL_ERROR) **/
  inline int set_function(UCL_Program &program, const char *function) {
    CUresult err=cuModuleGetFunction(&_kernel,program._module,function);
    if (err!=CUDA_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not find function: " << function
                << " in program.\n";
      UCL_GERYON_EXIT;
      #endif
      return UCL_FUNCTION_NOT_FOUND;
    }
    _cq=program._cq;
    return UCL_SUCCESS;
  }

  /// Set the kernel argument.
  /** If not a device pointer, this must be repeated each time the argument
    * changes
    * \note To set kernel parameter i (i>0), parameter i-1 must be set **/
  template <class dtype>
  inline void set_arg(const unsigned index, const dtype * const arg) {
    if (index==_num_args)
      add_arg(arg);
    else if (index<_num_args)
      #if CUDA_VERSION >= 4000
      _kernel_args[index]=arg;
      #else
      CU_SAFE_CALL(cuParamSetv(_kernel, _offsets[index], arg, sizeof(dtype)));
      #endif
    else
      assert(0==1); // Must add kernel parameters in sequential order
  }

  /// Set a geryon container as a kernel argument.
  template <class numtyp>
  inline void set_arg(const UCL_D_Vec<numtyp> * const arg)
    { set_arg(&arg->begin()); }

  /// Set a geryon container as a kernel argument.
  template <class numtyp>
  inline void set_arg(const UCL_D_Mat<numtyp> * const arg)
    { set_arg(&arg->begin()); }

  /// Set a geryon container as a kernel argument.
  template <class hosttype, class devtype>
  inline void set_arg(const UCL_Vector<hosttype, devtype> * const arg)
    { set_arg(&arg->device.begin()); }

  /// Set a geryon container as a kernel argument.
  template <class hosttype, class devtype>
  inline void set_arg(const UCL_Matrix<hosttype, devtype> * const arg)
    { set_arg(&arg->device.begin()); }

  /// Add a kernel argument.
  inline void add_arg(const CUdeviceptr* const arg) {
    #if CUDA_VERSION >= 4000
    _kernel_args[_num_args]=(void *)arg;
    #else
    void* ptr = (void*)(size_t)(*arg);
    _param_size = (_param_size + __alignof(ptr) - 1) & ~(__alignof(ptr) - 1);
    CU_SAFE_CALL(cuParamSetv(_kernel, _param_size, &ptr, sizeof(ptr)));
    _offsets.push_back(_param_size);
    _param_size+=sizeof(ptr);
    #endif
    _num_args++;
    if (_num_args>UCL_MAX_KERNEL_ARGS) assert(0==1);
  }

  /// Add a kernel argument.
  template <class dtype>
  inline void add_arg(const dtype* const arg) {
    #if CUDA_VERSION >= 4000
    _kernel_args[_num_args]=const_cast<dtype *>(arg);
    #else
    _param_size = (_param_size+__alignof(dtype)-1) & ~(__alignof(dtype)-1);
    CU_SAFE_CALL(cuParamSetv(_kernel,_param_size,(void*)arg,sizeof(dtype)));
    _offsets.push_back(_param_size);
    _param_size+=sizeof(dtype);
    #endif
    _num_args++;
    if (_num_args>UCL_MAX_KERNEL_ARGS) assert(0==1);
  }

  /// Add a geryon container as a kernel argument.
  template <class numtyp>
  inline void add_arg(const UCL_D_Vec<numtyp> * const arg)
    { add_arg(&arg->begin()); }

  /// Add a geryon container as a kernel argument.
  template <class numtyp>
  inline void add_arg(const UCL_D_Mat<numtyp> * const arg)
    { add_arg(&arg->begin()); }

  /// Add a geryon container as a kernel argument.
  template <class hosttype, class devtype>
  inline void add_arg(const UCL_Vector<hosttype, devtype> * const arg)
    { add_arg(&arg->device.begin()); }

  /// Add a geryon container as a kernel argument.
  template <class hosttype, class devtype>
  inline void add_arg(const UCL_Matrix<hosttype, devtype> * const arg)
    { add_arg(&arg->device.begin()); }

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue is used for the kernel execution **/
  inline void set_size(const size_t num_blocks, const size_t block_size) {
    _dimensions=1;
    _num_blocks[0]=num_blocks;
    _num_blocks[1]=1;
    _num_blocks[2]=1;
    #if CUDA_VERSION >= 4000
    _block_size[0]=block_size;
    _block_size[1]=1;
    _block_size[2]=1;
    #else
    CU_SAFE_CALL(cuFuncSetBlockShape(_kernel,block_size,1,1));
    #endif
  }

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue for the kernel is changed to cq **/
  inline void set_size(const size_t num_blocks, const size_t block_size,
                       command_queue &cq)
    { _cq=cq; set_size(num_blocks,block_size); }

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue is used for the kernel execution **/
  inline void set_size(const size_t num_blocks_x, const size_t num_blocks_y,
                       const size_t block_size_x, const size_t block_size_y) {
    _dimensions=2;
    _num_blocks[0]=num_blocks_x;
    _num_blocks[1]=num_blocks_y;
    _num_blocks[2]=1;
    #if CUDA_VERSION >= 4000
    _block_size[0]=block_size_x;
    _block_size[1]=block_size_y;
    _block_size[2]=1;
    #else
    CU_SAFE_CALL(cuFuncSetBlockShape(_kernel,block_size_x,block_size_y,1));
    #endif
  }

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue for the kernel is changed to cq **/
  inline void set_size(const size_t num_blocks_x, const size_t num_blocks_y,
                       const size_t block_size_x, const size_t block_size_y,
                       command_queue &cq)
    {_cq=cq; set_size(num_blocks_x, num_blocks_y, block_size_x, block_size_y);}

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue is used for the kernel execution **/
  inline void set_size(const size_t num_blocks_x, const size_t num_blocks_y,
                       const size_t block_size_x,
                       const size_t block_size_y, const size_t block_size_z) {
    _dimensions=2;
    _num_blocks[0]=num_blocks_x;
    _num_blocks[1]=num_blocks_y;
    _num_blocks[2]=1;
    #if CUDA_VERSION >= 4000
    _block_size[0]=block_size_x;
    _block_size[1]=block_size_y;
    _block_size[2]=block_size_z;
    #else
    CU_SAFE_CALL(cuFuncSetBlockShape(_kernel,block_size_x,block_size_y,
                                     block_size_z));
    #endif
  }

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue is used for the kernel execution **/
  inline void set_size(const size_t num_blocks_x, const size_t num_blocks_y,
                       const size_t block_size_x, const size_t block_size_y,
                       const size_t block_size_z, command_queue &cq) {
    _cq=cq;
    set_size(num_blocks_x, num_blocks_y, block_size_x, block_size_y,
             block_size_z);
  }

  /// Run the kernel in the default command queue
  inline void run() {
    #if CUDA_VERSION >= 4000
    CU_SAFE_CALL(cuLaunchKernel(_kernel,_num_blocks[0],_num_blocks[1],
                                _num_blocks[2],_block_size[0],_block_size[1],
                                _block_size[2],0,_cq,_kernel_args,nullptr));
    #else
    CU_SAFE_CALL(cuParamSetSize(_kernel,_param_size));
    CU_SAFE_CALL(cuLaunchGridAsync(_kernel,_num_blocks[0],_num_blocks[1],_cq));
    #endif
  }

  /// Clear any arguments associated with the kernel
  inline void clear_args() {
    _num_args=0;
    #if CUDA_VERSION < 4000
    _offsets.clear();
    _param_size=0;
    #endif
  }

  /// Return the default command queue/stream associated with this data
  inline command_queue & cq() { return _cq; }
  /// Change the default command queue associated with matrix
  inline void cq(command_queue &cq_in) { _cq=cq_in; }
  #include "ucl_arg_kludge.h"

 private:
  CUfunction _kernel;
  CUstream _cq;
  unsigned _dimensions;
  unsigned _num_blocks[3];
  unsigned _num_args;
  friend class UCL_Texture;

  #if CUDA_VERSION >= 4000
  unsigned _block_size[3];
  void * _kernel_args[UCL_MAX_KERNEL_ARGS];
  #else
  std::vector<unsigned> _offsets;
  unsigned _param_size;
  #endif
};

} // namespace

#endif

