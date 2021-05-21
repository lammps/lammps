/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author     : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project    : ivutils
 *
 *****************************************************************************/
/*s****************************************************************************
 * $Log: refobj.h,v $
 * Revision 1.2  2011/06/11 16:53:55  valuev
 * sync with LAMMPS
 *
 * Revision 1.1  2011/06/10 17:15:07  morozov
 * First Windows project with the correct directory structure
 *
 * Revision 1.15  2010/10/07 11:20:31  valuev
 * preliminary program restart
 *
 * Revision 1.14  2009/07/24 05:08:46  valuev
 * Sync with FDTD, added molecule setup
 *
 * Revision 1.33  2009/05/19 21:50:17  valuev
 * Added TestRay for plane
 *
 * Revision 1.32  2009/03/23 22:00:48  lesha
 * const mngptr &operator=(const mngarg<T> &arg) is added
 *
 * Revision 1.31  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.30  2009/01/21 09:28:15  lesha
 * refvector::clear is added
 *
 * Revision 1.29  2009/01/15 07:31:07  lesha
 * *** empty log message ***
 *
 * Revision 1.28  2009/01/14 10:02:36  lesha
 * operator [] is added to mngptr
 *
 * Revision 1.27  2008/04/29 01:23:59  lesha
 * nothing important
 *
 * Revision 1.26  2008/02/28 08:57:27  lesha
 * shptr::free() is made public
 *
 * Revision 1.25  2008/02/27 13:37:23  lesha
 * shptr is added
 *
 * Revision 1.24  2008/01/22 10:14:05  lesha
 * mngarg is added
 *
 * Revision 1.23  2007/08/08 10:55:37  lesha
 * constructor in refvector is correted
 *
 * Revision 1.22  2007/07/10 19:52:44  lesha
 * make gcc compilable
 *
 * Revision 1.21  2007/07/06 12:23:37  valuev
 * made compilable with icc 9
 *
 * Revision 1.20  2007/06/22 09:42:25  valuev
 * *** empty log message ***
 *
 * Revision 1.19  2007/06/17 00:51:44  lesha
 * refobj, refcounter :: reset is modified for ptr==this->ptr case
 *
 * Revision 1.18  2007/06/05 16:30:53  lesha
 * make gcc compilable
 *
 * Revision 1.17  2007/06/05 11:37:53  lesha
 * make gcc compilable
 *
 * Revision 1.16  2007/06/05 11:07:04  lesha
 * make gcc compilable
 *
 * Revision 1.15  2007/06/04 14:03:55  lesha
 * *** empty log message ***
 *
 * Revision 1.14  2007/05/31 18:00:42  lesha
 * *** empty log message ***
 *
 * Revision 1.13  2007/05/31 16:57:30  lesha
 * ref_sequence is added
 *
 * Revision 1.12  2007/05/31 01:25:01  lesha
 * new version of mng_ptr, pencil etc
 *
 * Revision 1.11  2007/02/20 10:26:11  valuev
 * added newlines at end of file
 *
 * Revision 1.10  2007/02/16 10:16:51  valuev
 * allowed array mngptr
 *
 * Revision 1.4  2007/02/16 09:40:32  valuev
 * Added Nudged Elastic Band saddle point search
 *
 * Revision 1.3  2006/12/20 14:29:33  valuev
 * Updated workflow, sync with FDTD
 *
 * Revision 1.9  2006/12/14 08:42:36  valuev
 * reformulated detector
 *  projectors, corrected open file limit control, tested Fourier sceleton
 *
 * Revision 1.8  2006/11/29 18:05:05  valuev
 * made the code compilable with g++
 *
 * Revision 1.7  2006/11/29 17:17:01  valuev
 * added using base_t::member for ANSI-compatibility
 *
 * Revision 1.6  2006/11/29 17:11:59  valuev
 * added includes
 *
 * Revision 1.5  2006/11/28 09:16:59  valuev
 * Fixed vectors storing managed pointers
 *
 * Revision 1.4  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
*******************************************************************************/
#ifndef _REFOBJ_H
#define _REFOBJ_H

# include <utility>
# include <vector>
# include <map>

using namespace std;

template<class T>
class mngarg: public pair<T *, int>{
public:
  typedef pair<T *, int> base_t;
  using base_t::second;
  using base_t::first;

  mngarg(T *ptr, int managed=0): pair<T*,int>(ptr,managed){}
  template<class A>
  mngarg(const mngarg<A> &arg): pair<T*,int>(arg.first,arg.second){}
};

template<class T>
mngarg<T> make_mngarg(T *ptr, int managed=1){
  return mngarg<T>(ptr,managed);
}

/// managed pointer
/// managed==0 do not delete
/// managed==1 delete
/// (NOT IMPLEMENTED) managed==2 copy and delete, requires copy constructor
/// flag 0x8 -- delete as array
template<class T>
class mngptr: public pair<T *, int>{
public:
  typedef pair<T *, int> base_t;
  typedef T *pointer;
  using base_t::second;
  using base_t::first;

  mngptr(T* ptr=nullptr, int managed=0): pair<T*,int>(ptr,managed){
    //if(managed==2)ptr= new T(*ptr);
  }
  mngptr(const mngarg<T> &arg): pair<T*,int>(arg.first,arg.second){}
  const mngptr &operator=(const mngarg<T> &arg){
    reset(arg.first,arg.second);
    return *this;
  }
  void reset(T* ptr=nullptr, int managed=0){
    if(second && first && first!=ptr){
      if(second&0x8)delete [] first;
      else delete first;
    }
    first=ptr;
    second=managed;
  }
  void reset(const mngarg<T> &arg){
    reset(arg.first,arg.second);
  }
  T* ptr() const {
    return first;
  }
  T* operator->() const {
    return first;
  }
  T& operator*() const{
    return *first;
  }
  T& operator[] (int i) const{
    return *(first+i);
  }
  int managed() const {
    return second;
  }
  ~mngptr(){
    reset();
  }
};


# if 0
template <template<class _type> class cont_tt, class T>
class refcontainer: public cont_tt< T * >{
protected:
  int man;
public:
  typedef cont_tt< T * > base_t;
  typedef typename base_t::iterator iterator;
  typedef typename base_t::const_iterator const_iterator;

  refcontainer(int smanaged=0):man(smanaged){}
  refcontainer(size_t n, int smanaged=0):base_t(n),man(smanaged){}

  void set_managed(int sman){
    man=sman;
  }

  ~refcontainer(){
    if(man){
      size_t i, n=base_t::size();
      for(i=0;i<n;i++)
        if((*this)[i]) delete (*this)[i];
    }
  }
};
# endif

template < class T >
class refvector: public std::vector< T * >{
protected:
  int man;
public:
  typedef vector<T*> base_t;
  typedef typename base_t::iterator iterator;
  typedef typename base_t::const_iterator const_iterator;

  refvector(int smanaged=0):man(smanaged){}
//  refvector(size_t n, int smanaged=0):base_t(n),man(smanaged){} // ambigious constructors
  refvector(size_t n, int smanaged):base_t(n),man(smanaged){}

  void set_managed(int sman){
    man=sman;
  }

  ~refvector(){
    clear();
  }

  void clear() {
    if(man){
      iterator it=base_t::begin();
      for(;it!=base_t::end();++it)
        if(*it)
          delete (*it);
    }
    base_t::clear();
  }

  iterator erase(iterator it){
    if(man && *it)
      delete (*it);
    return base_t::erase(it);
  }

};


template <class key_tt, class T >
class refmap: public std::map<key_tt, T * >{
protected:
  int man;
public:
  typedef std::map<key_tt, T * > base_t;
  typedef typename base_t::iterator iterator;
  typedef typename base_t::const_iterator const_iterator;

  refmap(int smanaged=0):man(smanaged){}
  refmap(size_t n, int smanaged=0):base_t(n),man(smanaged){}

  void set_managed(int sman){
    man=sman;
  }

  ~refmap(){
    clear();
  }

  void clear() {
    if(man){
      for(typename base_t::iterator i=base_t::begin();i!=base_t::end();++i)
        if(i->second) delete i->second;
    }
    base_t::clear();
  }

  iterator erase(iterator it){
    if(man && it->second)
      delete it->second;
    return base_t::erase(it);
  }
};


template<class T>
class delete_ptr{
public:
  void operator()(T *ptr){
    delete ptr;
  }
};

template<class T, class delete_t=delete_ptr<T> >
class shptr{
  template<class Y, class Z> friend class shptr;
  T *p;
  int *num; //if num==nullptr than p is not managed (as in mngptr)

  void set(T *p_, int managed){
    p=p_;
    if(p&&managed){
      num=new int;
      *num=1;
    }
    else num=nullptr;
  }
  template<class Y>
  void set(const Y &other){
    p=other.p;
    if(p){
      num=other.num;
      if(num)(*num)++;
    }
    else num=nullptr;
  }
  void set(const shptr &other){
    p=other.p;
    if(p){
      num=other.num;
      if(num)(*num)++;
    }
    else num=nullptr;
  }

public:
  shptr(T* p=nullptr, int managed=1){
    set(p,managed);
  }
  shptr(const mngarg<T> &arg){
    set(arg.first,arg.second);
  }
  template<class Y>
  shptr(const Y &other){
    set(other);
  }
  shptr(const shptr &other){
    set(other);
  }

  void reset(T *p_, int managed=1) {
    if(p!=p_){
      free();
      set(p_,managed);
    }
  }
  void reset(const shptr &other) {
    if(this!=&other){
      free();
      set(other);
    }
  }

  const shptr &operator=(T *p) {
    reset(p,0);
    return *this;
  }
  const shptr &operator=(const mngarg<T> &arg) {
    reset(arg.first,arg.second);
    return *this;
  }
  template<class Y>
  const shptr &operator=(const Y &other){
    reset(other);
    return *this;
  }
  const shptr &operator=(const shptr &other) {
    reset(other);
    return *this;
  }

  virtual ~shptr(){
    free();
  }
  void free(){
    if(p){
      if(num){
        (*num)--;
        if((*num)==0){
          delete_t()(p);
          delete num;
        }
        num=nullptr;
      }
      p=nullptr;
    }
  }

  bool valid() const {
    return p!=nullptr;
  }

  T* ptr() const {
    return p;
  }
  T* operator->() const {
    return p;
  }
  T& operator*() const{
    return *p;
  }
};




/*
class RefObject{
  void *ref_data;
  int  ref_count;
public:

protected:
  virtual void delete_data(void *data);
  virtual void *new_data();
  virtual void *copy_data(void *data);
}



class RefA: public RefObject{

public:
  refA(){
    ref_data = new A;
  }
  refA(const refA &other){
    Ref(other);
  }
  refA& operator=(const refA &other){
    if(ref_data != other.ref_data){
      Ref(other);
    }
    return *this;
  }
private:
  void delete_data(void *data){
    delete (A *)data;
  }
  void *new_data(){
    return (void *)new A;
  }
  void *copy_data(void *data){
    return (void *)(new A(*((A*)data)));
  }



}*/


#endif
