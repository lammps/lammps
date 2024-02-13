/***************************************************************************
                              ucl_s_obj_help.h
                             -------------------
                               W. Michael Brown

  Helper routines for allocating memory for s-objects and performing
  host/device updates. (Different routines depending on whether the
  same type is used on the host and device).

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Mon May 14 2012
    copyright            : (C) 2012 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
    This software is distributed under the Simplified BSD License.
   ----------------------------------------------------------------------- */

template <int st> struct _ucl_s_obj_help;

// Host and device containers are same type
// -- Don't need casting buffers
// -- Can potentially use same memory if shared by accelerator
template <> struct _ucl_s_obj_help<1> {
  template <class t1, class t2, class t3>
    static inline int alloc(t1 &host, t2 &device, t3 & /*_buffer*/,
                          const int cols, UCL_Device &acc,
                          const enum UCL_MEMOPT kind1,
                          const enum UCL_MEMOPT kind2) {
    int e1;
    if (acc.shared_memory()) {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 1S\n";
      #endif
      e1=host.alloc(cols,acc,kind1,kind2);
      if (e1!=UCL_SUCCESS)
        return e1;
      device.view(host);
      return UCL_SUCCESS;
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 1NS\n";
      #endif
      e1=host.alloc(cols,acc,kind1);
      if (e1!=UCL_SUCCESS)
        return e1;
      return device.alloc(cols,acc,kind2);
    }
  }

  template <class t1, class t2, class t3, class mat_type>
  static inline int alloc(t1 &host, t2 &device, t3 &/*_buffer*/,
                          const int cols, mat_type &cq,
                          const enum UCL_MEMOPT kind1,
                          const enum UCL_MEMOPT kind2) {
    int e1;
    if (cq.shared_mem_device()) {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 2S\n";
      #endif
      e1=host.alloc(cols,cq,kind1,kind2);
      if (e1!=UCL_SUCCESS)
        return e1;
      device.view(host);
      return UCL_SUCCESS;
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 2NS\n";
      #endif
      e1=host.alloc(cols,cq,kind1);
      if (e1!=UCL_SUCCESS)
        return e1;
      return device.alloc(cols,cq,kind2);
    }
  }

  template <class t1, class t2, class t3>
  static inline int alloc(t1 &host, t2 &device, t3 &/*_buffer*/,
                          const int rows, const int cols, UCL_Device &acc,
                          const enum UCL_MEMOPT kind1,
                          const enum UCL_MEMOPT kind2) {
    int e1;
    if (acc.shared_memory()) {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 3S\n";
      #endif
      e1=host.alloc(rows,cols,acc,kind1,kind2);
      if (e1!=UCL_SUCCESS)
        return e1;
      device.view(host);
      return UCL_SUCCESS;
    } else {
      e1=host.alloc(rows,cols,acc,kind1);
      if (e1!=UCL_SUCCESS)
        return e1;
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 3NS\n";
      #endif
      return device.alloc(rows,cols,acc,kind2);
    }
  }

  template <class t1, class t2, class t3, class mat_type>
  static inline int alloc(t1 &host, t2 &device, t3 &/*_buffer*/,
                          const int rows, const int cols, mat_type &cq,
                          const enum UCL_MEMOPT kind1,
                          const enum UCL_MEMOPT kind2) {
    int e1;
    if (cq.shared_mem_device()) {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 4S\n";
      #endif
      e1=host.alloc(rows,cols,cq,kind1,kind2);
      if (e1!=UCL_SUCCESS)
        return e1;
      device.view(host);
      return UCL_SUCCESS;
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 4NS\n";
      #endif
      e1=host.alloc(rows,cols,cq,kind1);
      if (e1!=UCL_SUCCESS)
        return e1;
      return device.alloc(rows,cols,cq,kind2);
    }
  }

  template <class t1, class t2, class t3>
    static inline void copy(t1 &dst, t2 &src, t3 & /*buffer*/, const bool async) {
    ucl_copy(dst,src,async);
  }

  template <class t1, class t2, class t3>
    static inline void copy(t1 &dst, t2 &src, t3 & /*buffer*/, command_queue &cq) {
    ucl_copy(dst,src,cq);
  }

  template <class t1, class t2, class t3>
    static inline void copy(t1 &dst, t2 &src, const int cols, t3 & /*buffer*/, const bool async) {
    ucl_copy(dst,src,cols,async);
  }

  template <class t1, class t2, class t3>
    static inline void copy(t1 &dst, t2 &src, const int cols, t3 & /*buffer*/, command_queue &cq) {
    ucl_copy(dst,src,cols,cq);
  }

  template <class t1, class t2, class t3>
    static inline void copy(t1 &dst, t2 &src, const int rows, const int cols, t3 & /*buffer*/, const bool async) {
    ucl_copy(dst,src,rows,cols,async);
  }

  template <class t1, class t2, class t3>
  static inline void copy(t1 &dst, t2 &src, const int rows, const int cols, t3 & /*buffer*/, command_queue &cq) {
    ucl_copy(dst,src,rows,cols,cq);
  }

  template <class t1, class t2, class t3>
    static inline int dev_resize(t1 &device, t2 &host, t3 & /*buff*/,const int cols) {
    if (device.kind()==UCL_VIEW) {
      device.view(host);
      return UCL_SUCCESS;
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 5S\n";
      #endif
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 5NS\n";
      #endif
      return device.resize(cols);
    }
  }

  template <class t1, class t2, class t3>
  static inline int dev_resize(t1 &device, t2 &host, t3 &/*buff*/, const int rows,
                               const int cols) {
    if (device.kind()==UCL_VIEW) {
      device.view(host);
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 6S\n";
      #endif
      return UCL_SUCCESS;
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 6NS\n";
      #endif
      return device.resize(rows,cols);
    }
  }
};

// Host and device containers are different types
template <int st> struct _ucl_s_obj_help {
  template <class t1, class t2, class t3>
  static inline int alloc(t1 &host, t2 &device, t3 &_buffer,
                          const int cols, UCL_Device &acc,
                          const enum UCL_MEMOPT kind1,
                          const enum UCL_MEMOPT kind2) {
    int e1;
    e1=host.alloc(cols,acc,UCL_NOT_PINNED);
    if (e1!=UCL_SUCCESS)
      return e1;

    if (acc.shared_memory()) {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 7S\n";
      #endif
      e1=_buffer.alloc(cols,acc,kind1,kind2);
      if (e1!=UCL_SUCCESS)
        return e1;
      device.view(_buffer);
      return UCL_SUCCESS;
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 7NS\n";
      #endif
      e1=_buffer.alloc(cols,acc,kind1);
      if (e1!=UCL_SUCCESS)
        return e1;
      return device.alloc(cols,acc,kind2);
    }
  }

  template <class t1, class t2, class t3, class mat_type>
  static inline int alloc(t1 &host, t2 &device, t3 &_buffer,
                          const int cols, mat_type &cq,
                          const enum UCL_MEMOPT kind1,
                          const enum UCL_MEMOPT kind2) {
    int e1;
    e1=host.alloc(cols,cq,UCL_NOT_PINNED);
    if (e1!=UCL_SUCCESS)
      return e1;
    if (cq.shared_mem_device()) {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 8S\n";
      #endif
      e1=_buffer.alloc(cols,cq,kind1,kind2);
      if (e1!=UCL_SUCCESS)
        return e1;
      device.view(_buffer);
      return UCL_SUCCESS;
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 8NS\n";
      #endif
      e1=_buffer.alloc(cols,cq,kind1);
      if (e1!=UCL_SUCCESS)
        return e1;
      return device.alloc(cols,cq,kind2);
    }
  }

  template <class t1, class t2, class t3>
  static inline int alloc(t1 &host, t2 &device, t3 &_buffer,
                          const int rows, const int cols, UCL_Device &acc,
                          const enum UCL_MEMOPT kind1,
                          const enum UCL_MEMOPT kind2) {
    int e1;
    e1=host.alloc(rows,cols,acc,UCL_NOT_PINNED);
    if (e1!=UCL_SUCCESS)
      return e1;

    if (acc.shared_memory()) {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 9S\n";
      #endif
      e1=_buffer.alloc(rows,cols,acc,kind1,kind2);
      if (e1!=UCL_SUCCESS)
        return e1;
      device.view(_buffer);
      return UCL_SUCCESS;
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 9NS\n";
      #endif
      e1=_buffer.alloc(rows,cols,acc,kind1);
      if (e1!=UCL_SUCCESS)
        return e1;
      return device.alloc(rows,cols,acc,kind2);
    }
  }

  template <class t1, class t2, class t3, class mat_type>
  static inline int alloc(t1 &host, t2 &device, t3 &_buffer,
                          const int rows, const int cols, mat_type &cq,
                          const enum UCL_MEMOPT kind1,
                          const enum UCL_MEMOPT kind2) {
    int e1;
    e1=host.alloc(rows,cols,cq,UCL_NOT_PINNED);
    if (e1!=UCL_SUCCESS)
      return e1;
    if (cq.shared_mem_device()) {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 10S\n";
      #endif
      e1=_buffer.alloc(rows,cols,cq,kind1,kind2);
      if (e1!=UCL_SUCCESS)
        return e1;
      device.view(_buffer);
      return UCL_SUCCESS;
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 10NS\n";
      #endif
      e1=_buffer.alloc(rows,cols,cq,kind1);
      if (e1!=UCL_SUCCESS)
        return e1;
      return device.alloc(rows,cols,cq,kind2);
    }
  }

  template <class t1, class t2, class t3>
  static inline void copy(t1 &dst, t2 &src, t3 &buffer, const bool async) {
    ucl_cast_copy(dst,src,buffer,async);
  }

  template <class t1, class t2, class t3>
  static inline void copy(t1 &dst, t2 &src, t3 &buffer, command_queue &cq) {
    ucl_cast_copy(dst,src,buffer,cq);
  }

  template <class t1, class t2, class t3>
  static inline void copy(t1 &dst, t2 &src, const int cols, t3 &buffer,
                          const bool async) {
    ucl_cast_copy(dst,src,cols,buffer,async);
  }

  template <class t1, class t2, class t3>
  static inline void copy(t1 &dst, t2 &src, const int cols, t3 &buffer,
                          command_queue &cq) {
    ucl_cast_copy(dst,src,cols,buffer,cq);
  }

  template <class t1, class t2, class t3>
  static inline void copy(t1 &dst, t2 &src, const int rows, const int cols,
                          t3 &buffer, const bool async) {
    ucl_cast_copy(dst,src,rows,cols,buffer,async);
  }

  template <class t1, class t2, class t3>
  static inline void copy(t1 &dst, t2 &src, const int rows, const int cols,
                          t3 &buffer, command_queue &cq) {
    ucl_cast_copy(dst,src,rows,cols,buffer,cq);
  }

  template <class t1, class t2, class t3>
  static inline int dev_resize(t1 &device, t2 & /*host*/, t3 &buff,const int cols) {
    int err=buff.resize(cols);
    if (err!=UCL_SUCCESS)
      return err;

    if (device.kind()==UCL_VIEW) {
      device.view(buff);
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 11S\n";
      #endif
      return UCL_SUCCESS;
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 11NS\n";
      #endif
      return device.resize(cols);
    }
  }

  template <class t1, class t2, class t3>
  static inline int dev_resize(t1 &device, t2 &/*host*/, t3 &buff, const int rows,
                               const int cols) {
    int err=buff.resize(rows,cols);
    if (err!=UCL_SUCCESS)
      return err;

    if (device.kind()==UCL_VIEW) {
      device.view(buff);
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 12S\n";
      #endif
      return UCL_SUCCESS;
    } else {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_ALLOC 12NS\n";
      #endif
      return device.resize(rows,cols);
    }
  }

};

