/*e***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author  : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project  : ivutils
 *
 *
 *****************************************************************************/
/*e****************************************************************************
 * $Log: pairhash.h,v $
 * Revision 1.3  2011/06/11 18:18:50  morozov
 * USER-AWPMD compiles on Linux now!
 *
 * Revision 1.2  2011/06/11 16:53:55  valuev
 * sync with LAMMPS
 *
 * Revision 1.1  2011/06/10 17:15:07  morozov
 * First Windows project with the correct directory structure
 *
 * Revision 1.27  2011/06/09 22:55:08  valuev
 * norm matrices
 *
 * Revision 1.26  2011/06/07 17:43:00  valuev
 * added Y derivatives
 *
 * Revision 1.25  2011/05/24 19:54:32  valuev
 * fixed sqmatrix::iterator
 *
 * Revision 1.24  2011/05/21 23:06:49  valuev
 * Norm matrix transform to pysical variables
 *
 * Revision 1.23  2009/09/24 10:06:38  valuev
 * moved matrix printing to template function, reproducing old TB calculations
 *
 * Revision 1.22  2009/02/10 14:20:45  valuev
 * sync with FDTD project
 *
 * Revision 1.4  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.21  2008/08/27 13:34:32  valuev
 * made icc-compilable
 *
 * Revision 1.20  2008/08/25 21:06:11  valuev
 * moved using delaration to public
 *
 * Revision 1.19  2008/07/23 16:55:05  morozov
 * *** empty log message ***
 *
 * Revision 1.18  2008/07/23 16:21:52  morozov
 * Corrected Makefile for unix compilation of tcpengine
 *
 * Revision 1.17  2008/07/18 18:15:31  morozov
 * *** empty log message ***
 *
 * Revision 1.16  2008/07/02 13:11:32  valuev
 * new C60+O2 experiments
 *
 * Revision 1.15  2008/06/24 08:50:00  valuev
 * made icc-compilable
 *
 * Revision 1.14  2008/06/24 08:39:57  valuev
 * added ESSL support to TB
 *
 * Revision 1.13  2008/05/29 14:47:33  valuev
 * made icc-compilable
 *
 * Revision 1.12  2008/05/14 17:17:22  morozov
 * Passed 2- and 3-electron test. Added Norm matrix.
 *
 * Revision 1.11  2008/05/05 17:27:43  morozov
 * cvector_3.h is the new header for class cVector_3. Old one is moved to cvector_3old.h
 * class hmatrix is added to pairhash.h
 * wavepackets.h contains class WavePacket
 *
 * Revision 1.10  2008/04/21 22:42:30  valuev
 * *** empty log message ***
 *
 * Revision 1.9  2008/04/15 13:11:41  valuev
 * Added antisymmetrized wave packets
 *
 * Revision 1.8  2008/02/28 13:26:04  valuev
 * vasp scanner
 *
 * Revision 1.7  2007/12/13 19:48:59  valuev
 * added newlines
 *
 * Revision 1.6  2006/12/20 14:29:33  valuev
 * Updated workflow, sync with FDTD
 *
 * Revision 1.3  2006/10/27 20:41:01  valuev
 * Added detectors sceleton. Updated some of ivutils from MD project.
 *
 * Revision 1.5  2006/09/26 10:59:42  valuev
 * Added nonorthogonal TB (Menon-Subbaswamy)
 *
 * Revision 1.4  2006/07/21 16:22:03  valuev
 * Added Tight Binding for graphite+O
 *
 * Revision 1.3  2006/04/26 12:12:01  valuev
 * Fixed Neighbour Lists (double-single counting), added twostep NL scheme, added Step2 to mdtutorial (use of template potentials), added DelAtom to mdStructure
 *
 * Revision 1.2  2005/12/09 21:06:38  valuev
 * Added neighbour list to mdPotential interface.
 * Added missing files to ivutils directory.
 * Added mdtutorial and step1 project
 *
 * Revision 1.1  2005/12/02 18:51:06  valuev
 * added  HEAD project tree
 *
 * Revision 1.1  2005/11/30 23:36:11  valuev
 * put ivutils to cvs on biolab1.mipt.ru
 *
 * Revision 1.1  2005/11/30 23:15:43  valuev
 * put ivutils on cvs biolab1.mipt.ru
 *
 *
*******************************************************************************/
# ifndef PAIRHASH_H
# define PAIRHASH_H

/*e @file pairhash.h @brief pair hash table
*/ 


/*r @file pairhash.h @brief работа с хеш-таблицами парных величин 
*/ 

# include "refobj.h"


///\en Rectangular matrix
template <class T>
class recmatrix{
protected:
  mngptr<T> parr;
public:
  class iterator{
    friend class recmatrix<T>;
    T *ptr;
    size_t incr;
    iterator(const recmatrix<T> *parent,size_t first_,size_t second_, bool inc_first=false){
      ptr=parent->arr+parent->index(first_,second_);
      incr=inc_first ? parent->sizex : 1 ;
    };
    iterator(T *ptr_,size_t incr_):ptr(ptr_),incr(incr_){}
  public: 
    iterator(const iterator &other):ptr(other.ptr),incr(other.incr){
    }
    iterator():ptr(NULL),incr(0){}
    iterator &operator++(){ // prefix
      ptr+=incr;
      return *this;
    }
    iterator operator++(int){ // postfix
      iterator tmp=*this;
      ++*this;
      return tmp;
    }
    iterator operator+(int delta) const { 
      return iterator(ptr+delta*incr,incr);
    }
    bool operator!=(const iterator &other) const {
      if(ptr!=other.ptr)return true;
      else return false;
    }
    T &operator*() const {
      return *ptr;
    }
  };


  T *arr;
  size_t sizex, sizey;

  //e default constructor
  recmatrix(): parr(NULL,1) {
    sizey=sizex=0;
    arr=NULL;
  }

  //e copy constructor: makes a managed copy
  recmatrix(const recmatrix &other):sizex(0),sizey(0),arr(NULL){
    *this=other;
  }

  recmatrix &operator=(const recmatrix &other){
    if(this!=&other){
      if(other.sizex*other.sizey<=sizex*sizey)
        init(other.sizex,other.sizey,-1); // keeping old array
      else
        init(other.sizex,other.sizey,1);
      size_t n=get_datasize(sizex,sizey);
      for(size_t i=0;i<n;i++)
        arr[i]=other.arr[i];  
    }
    return *this;
  }

  virtual size_t get_datasize(size_t nx, size_t ny) const{
    return nx*ny;
  }

  // i is y (row number), j is x (column number)
  size_t index(size_t i, size_t j) const {
    return sizey*i+j;
  }

  T &operator()(size_t i,size_t j){
    return arr[index(i,j)];
  }

  T operator()(size_t i,size_t j) const {
    return arr[index(i,j)];
  }


  void set(long i,long j,const T& val){
    (*this)(i,j)=val;
  }


  virtual int init(size_t nx, size_t ny, int smanaged=-1){
    int managed=parr.managed();
    if(managed && (sizex!=nx || sizey!=ny)){  
      parr.reset(NULL,0);
    }
    if(smanaged>=0){ // for changing the managed flag?
      parr.reset(parr.ptr(),smanaged ? smanaged|0x8 : 0  );
      managed=smanaged;
    }
    if(sizex==nx && sizey==ny) // no need to allocate
      return 1;
    sizex=nx;
    sizey=ny;
    if(managed){
      if(sizex>0 && sizey>0)
        parr.reset(new T[get_datasize(sizex,sizey)],managed|0x8);
      arr=parr.ptr();
    }
    return 1;
  } 

  recmatrix(size_t nx, size_t ny):sizex(0), sizey(0){ 
    init(nx,ny,1);
  }

  //e initializes by unmanaged pointer
  recmatrix(size_t nx, size_t ny , T *ptr):parr(ptr,0),sizex(nx), sizey(ny) { 
    init(nx,ny);
  }

  //e attaches to new pointer and sets unmanaged size
  void AttachTo(size_t nx,size_t ny, T *ptr){
    init(0,0);
    sizex=nx;
    sizey=ny;
    parr.reset(ptr,0);
  }

  void Set(const T &val){
    size_t i, n=get_datasize(sizex,sizey);
    for(i=0;i<n;i++)arr[i]=val;
  }
  
  void SetDiag(const T &val){
    size_t i, size=(sizex>sizey? sizey: sizex);
    for(i=0;i<size;i++){
      (*this)(i,i)=val;
    }
  }
 
  /// returns iterator with fixed first index to iterate through matrix line elements
  iterator fix_first(size_t first, size_t second) const {
    return iterator(this,first,second,false);
  }

  /// returns iterator with fixed second index to iterate through matrix column elements
  iterator fix_second(size_t first, size_t second) const {
    return iterator(this,first,second, true);
  }
 
};


//e square matrix
template <class T>
class sqmatrix: public recmatrix<T> {
  
public:

  size_t size;

  //e default constructor
  sqmatrix(){}

  //e copy constructor: makes a managed copy
  sqmatrix(const sqmatrix &other):size(0){
    *this=other;
  }

  sqmatrix &operator=(const sqmatrix &other){
    if(this!=&other){
      *((recmatrix<T> *)this)=*((recmatrix<T> *)&other);
      size=other.size;
    }
    return *this;
  }

  virtual size_t get_datasize(size_t n) const{
    return n*n;
  }

 
  virtual int init(size_t n, int smanaged=-1){
    size=n;
    return recmatrix<T>::init(n,n,smanaged);
  } 

  sqmatrix(size_t n):size(0){ 
    init(n,1);
  }

  //e initializes by unmanaged pointer
  sqmatrix(size_t n, T * /* ptr */):size(n){ 
    init(n);
  }

  //e attaches to new pointer and sets unmanaged size
  void AttachTo(size_t n, T *ptr){
    init(0);
    size=n;
    recmatrix<T>::parr.reset(ptr,0);
  }

 
};



# if 0
//e square matrix
template <class T>
class sqmatrix{
  mngptr<T> parr;
public:
  class iterator{
    friend class sqmatrix<T>;
    T *ptr;
    size_t incr;
    iterator(const sqmatrix<T> *parent,size_t first_,size_t second_, bool inc_first=false){
      ptr=parent->arr+parent->index(first_,second_);
      incr=inc_first ? parent->size : 1 ;
    };
    iterator(T *ptr_,size_t incr_):ptr(ptr_),incr(incr_){}
  public: 
    iterator(const iterator &other):ptr(other.ptr),incr(other.incr){
    }
    iterator():ptr(NULL),incr(0){}
    iterator &operator++(){ // prefix
      ptr+=incr;
      return *this;
    }
    iterator operator++(int){ // postfix
      iterator tmp=*this;
      ++*this;
      return tmp;
    }
    iterator operator+(int delta) const { 
      return iterator(ptr+delta*incr,incr);
    }
    bool operator!=(const iterator &other) const {
      if(ptr!=other.ptr)return true;
      else return false;
    }
    T &operator*() const {
      return *ptr;
    }
  };


  T *arr;
  size_t size;

  //e default constructor
  sqmatrix(): parr(NULL,1) {
    size=0;
    arr=NULL;
  }

  //e copy constructor: makes a managed copy
  sqmatrix(const sqmatrix &other):size(0),arr(NULL){
    *this=other;
  }

  sqmatrix &operator=(const sqmatrix &other){
    if(this!=&other){
      if(other.size<=size)
        init(other.size,-1); // keeping old array
      else
        init(other.size,1);
      size_t n=get_datasize(size);
      for(size_t i=0;i<n;i++)
        arr[i]=other.arr[i];  
    }
    return *this;
  }

  virtual size_t get_datasize(size_t n) const{
    return n*n;
  }

  size_t index(size_t i, size_t j) const {
    return size*i+j;
  }

  T &operator()(size_t i,size_t j){
    return arr[index(i,j)];
  }

  T operator()(size_t i,size_t j) const {
    return arr[index(i,j)];
  }


  void set(long i,long j,const T& val){
    (*this)(i,j)=val;
  }


  virtual int init(size_t n, int smanaged=-1){
    int managed=parr.managed();
    if(managed && size!=n){  
      parr.reset(NULL,0);
    }
    if(smanaged>=0){ // for changing the managed flag?
      parr.reset(parr.ptr(),smanaged ? smanaged|0x8 : 0  );
      managed=smanaged;
    }
    if(size==n) // no need to allocate
      return 1;
    size=n;
    if(managed){
      if(size>0)
        parr.reset(new T[get_datasize(size)],managed|0x8);
      arr=parr.ptr();
    }
    return 1;
  } 

  sqmatrix(size_t n):size(0){ 
    init(n,1);
  }

  //e initializes by unmanaged pointer
  sqmatrix(size_t n, T *ptr):parr(ptr,0),size(n){ 
    init(n);
  }

  //e attaches to new pointer and sets unmanaged size
  void AttachTo(size_t n, T *ptr){
    init(0);
    size=n;
    parr.reset(ptr,0);
  }

  void Set(const T &val){
    size_t i, n=get_datasize(size);
    for(i=0;i<n;i++)arr[i]=val;
  }
  
  void SetDiag(const T &val){
    size_t i;
    for(i=0;i<size;i++){
      (*this)(i,i)=val;
    }
  }
 
  /// returns iterator with fixed first index to iterate through matrix line elements
  iterator fix_first(size_t first, size_t second) const {
    return iterator(this,first,second,false);
  }

  /// returns iterator with fixed second index to iterate through matrix column elements
  iterator fix_second(size_t first, size_t second) const {
    return iterator(this,first,second, true);
  }
 
};

# endif 

 //e prints the matrix into a file
template< class matrix_t>
int fileout(FILE *f, const matrix_t &matr, const char *elm_fmt, const char *elm_sep=" ", const char *line_sep="\n"){
  size_t i, j; 
  int res=0;
  for(i=0;i<matr.size;i++){
    for(j=0;j<matr.size;j++){
      res+=fprintf(f,elm_fmt,matr(i,j));
      fprintf(f,elm_sep);
    }
    fprintf(f,line_sep);
  }
  return res;
}


//e symmetric matrix
template <class T>
class smatrix: public sqmatrix<T>{
  typedef  sqmatrix<T> base_t;
public:
 
  virtual size_t get_datasize(size_t n) const{
    return n*(n+1)/2;
  }
  
  size_t index(size_t i, size_t j) const {
    if(i>=j)
      return (2*base_t::size-j-1)*j/2+i;
    else 
      return (2*base_t::size-i-1)*i/2+j;
  }

  T &operator()(size_t i,size_t j){
    return base_t::arr[index(i,j)];
  }

  T operator()(size_t i,size_t j) const {
    return base_t::arr[index(i,j)];
  }

  void set(long i,long j,const T& val){
    (*this)(i,j)=val;
  }

  void SetDiag(const T &val){
    size_t i;
    for(i=0;i<base_t::size;i++){
      (*this)(i,i)=val;
    }
  }

  smatrix(){}

  //e copy constructor: makes a managed copy
  smatrix(const smatrix &other):sqmatrix<T>(other){}

  smatrix &operator=(const smatrix &other){
    return (smatrix&)( *(sqmatrix<T> *)this = (sqmatrix<T> &)other );
  }

  smatrix(size_t n):sqmatrix<T>(n){}

  //e initializes by unmanaged pointer
  smatrix(size_t n, T *ptr):sqmatrix<T>(n,ptr){}

};



//e Hermitian matrix
template<class T>
class hmatrix : public smatrix<T>
{
public:
  using smatrix<T>::arr;
  using smatrix<T>::size;
  hmatrix() : smatrix<T>() {}

  hmatrix(const smatrix<T> &other) : smatrix<T>(other) {}
  
  //e copy constructor: makes a managed copy
  hmatrix(const hmatrix &other): smatrix<T>(other){}

  hmatrix &operator=(const hmatrix &other) {
    return (hmatrix&)( *(smatrix<T>*)this = (smatrix<T>&)other );
  }

  hmatrix(size_t n) : smatrix<T>(n) {}

  hmatrix(size_t n, T *ptr) : smatrix<T>(n, ptr) {}

  T operator()(size_t i, size_t j) const {
    if(i<=j) return arr[(2*size-i-1)*i/2+j];
    else return conj( arr[(2*size-j-1)*j/2+i] );
  }

  void set(long i,long j,const T& val){
    if(i<=j) arr[(2*size-i-1)*i/2+j] = val;
    else arr[(2*size-j-1)*j/2+i] = conj(val);
  }

};


//e Basic pair hash class
template <class T>
class PairHash{
public:
  //e find the value with indexes i, j
  //e @return 0 if not found, 1 otherwise
  //e if retval is not NULL, puts the found value there
  virtual int Find(long i, long j, T *retval=NULL)=0;
  virtual int Find(long i, long j, T **retval=NULL)=0;
  virtual int Del(long i, long j)=0;
  virtual int Put(long i, long j, const T *value)=0;
  virtual int Put(long i, long j, const T& value)=0;
  virtual int Clear()=0;
  virtual ~PairHash(){}
};

//e Hash with symmetric matrix
template <class T>
class PairHashM: public PairHash<T>{
  smatrix<long> indm;
  T *arr;
  int as;
public:
  PairHashM(long n, int antisymmetric =0):indm(n),as(antisymmetric){
    indm.Set(-1);
    arr= new T[n*(n+1)/2];
  }
  int Find(long i, long j, T *retval=NULL){
    long ind=indm(i,j);
    if(ind>=0){
      if(retval){
        if(as && i<j)*retval=-arr[ind];
        else *retval=arr[ind];
      }
      return 1;
    }
    return 0;
  }
  int Find(long i, long j, T **retval){
    long ind=indm(i,j);
    if(ind>=0){
      *retval=&arr[ind];
      if(as && i<j)return -1;
      return 1;
    }
    return 0;
  }


  int Del(long i, long j){
    indm(i,j)=-1;
    return 1;
  }

  int Put(long i, long j, const T *value){
    long ind=indm.index(i,j);
    indm.arr[ind]=ind;
    arr[ind]=*value;
    if(as && i<j)arr[ind]=-arr[ind];
    return 1;
  }

  int Put(long i, long j, const T& value){
    long ind=indm.index(i,j);
    indm.arr[ind]=ind;
    arr[ind]=value;
    if(as && i<j)arr[ind]=-arr[ind];
    return 1;
  }

  int Clear(){
    indm.Set(-1);
    return 1;
  }
  
  virtual ~PairHashM(){
    delete [] arr;
  }
};



# endif
