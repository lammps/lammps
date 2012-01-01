/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *****************************************************************************/

/*s****************************************************************************
 * $Log: vector_3.h,v $
 * Revision 1.2  2011/06/11 16:53:55  valuev
 * sync with LAMMPS
 *
 * Revision 1.1  2011/06/10 17:15:07  morozov
 * First Windows project with the correct directory structure
 *
 * Revision 1.22  2010/11/17 02:13:32  valuev
 * Added analysis phase, fixed some memory leaks
 *
 * Revision 1.21  2009/09/09 14:43:35  valuev
 * fixed trajreader
 *
 * Revision 1.20  2009/07/24 05:08:46  valuev
 * Sync with FDTD, added molecule setup
 *
 * Revision 1.33  2009/06/01 13:01:42  valuev
 * Added ShearBox
 *
 * Revision 1.32  2009/05/28 07:49:00  valuev
 * updated chemosensor
 *
 * Revision 1.31  2009/05/22 08:53:52  valuev
 * new uiExperiment interface
 *
 * Revision 1.30  2009/02/22 10:49:56  lesha
 * Vector_Nt constructors are modified
 *
 * Revision 1.29  2009/02/04 14:23:31  valuev
 * fixed bug in maxabscoord and minabscoord functions
 *
 * Revision 1.28  2009/02/04 10:56:30  lesha
 * SINGLE_PRECISION is recovered
 *
 * Revision 1.27  2009/02/03 00:47:35  lesha
 * *** empty log message ***
 *
 * Revision 1.26  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.16  2008/08/18 21:40:09  valuev
 * added Gurski-Krasko potential
 *
 * Revision 1.15  2008/05/05 17:27:43  morozov
 * cvector_3.h is the new header for class cVector_3. Old one is moved to cvector_3old.h
 * class hmatrix is added to pairhash.h
 * wavepackets.h contains class WavePacket
 *
 * Revision 1.14  2008/04/22 12:44:17  valuev
 * made gcc 4.12 compilable
 *
 * Revision 1.13  2008/04/21 23:13:44  valuev
 * made gcc 4.12 compilable
 *
 * Revision 1.12  2008/04/21 22:42:30  valuev
 * *** empty log message ***
 *
 * Revision 1.11  2008/04/15 13:11:41  valuev
 * Added antisymmetrized wave packets
 *
 
 *
*******************************************************************************/
/*r @file vector_3.h @brief работа с трехмерными векторами 
*/ 

# ifndef VECTOR_3_H

# define VECTOR_3_H

# define _USE_MATH_DEFINES

# include <iostream>
# include <cmath>
# include <limits>
# include <stdlib.h>

// some compilers don't define PI!
# ifndef M_PI
# define M_PI 3.1415926535897932385
# endif

# ifndef fmod
//r деление по модулю чисел с плавающей точкой
# define fmod(a,b)  ((a)-((long)((a)/(b))*(b)))
# endif

using namespace std;

#ifndef SINGLE_PRECISION
typedef double vec_type;
#else
typedef float vec_type;
#endif

//e "infinitely" large number in Vector_3 sense
//r "бесконечно большое" число в смысле Vector_3
//# ifndef SINGLE_PRECISION
//# define VEC_INFTY 1.e20
//# else
//# define VEC_INFTY 1.e20f
//# endif
#define VEC_INFTY numeric_limits<vec_type>::max()

//e "infinitely" small number in Vector_3 sense (used for vector comparisons)
//r "бесконечно малое" число в смысле Vector_3 (используется для сравнений)
//# ifndef SINGLE_PRECISION
//# define VEC_ZERO 1.e-20
//# else
//# define VEC_ZERO 1.e-20f
//# endif
#define VEC_ZERO 512*numeric_limits<vec_type>::epsilon()

//r N-мерный вектор типа T, с некоторыми полезными операциями
template <class T, const int N=3> 
struct Vector_Nt {
  typedef T value_type;
  
  T v[N];

  Vector_Nt(const T &a=0){
    for (int i=0; i<N; i++)
      v[i]=a;
  }

  
  explicit Vector_Nt(const T &a1, const T &a2) {
    if(N>0)v[0]=a1;
    if(N>1)v[1]=a2;
    for(int i=2;i<N;i++)v[i]=0;
  }

  explicit Vector_Nt(const T &s1, const T &s2, const T& s3) {
    if(N>0)v[0]=s1;
    if(N>1)v[1]=s2;
    if(N>2)v[2]=s3;
    for(int i=3;i<N;i++)v[i]=0;
  }
  //e construct from input iterator (or array)
  template <class A>
  Vector_Nt(const A *beg) {
    for (int i=0; i<N; i++, ++beg)
      v[i]=*beg;
  };

  //e construct from another vector
  template <class A>
  Vector_Nt(const Vector_Nt<A,N>& v1) {
    for (int i=0; i<N; i++) v[i]=v1[i];
  };

  //r Копирует содержимое вектора в it
  template <class A>
  void copy_to(A *beg) const {  
    for (int i=0; i<N; i++, ++beg)
      *beg=v[i];
  };

  //r получение элемента 
  inline T& operator[](int i) const { return (T&)v[i]; };


  inline Vector_Nt& operator=(const Vector_Nt &vect){
    for (int i=0; i<N; i++)
      v[i]=vect.v[i];
    return *this;
  }

  //r присваивает всем компонентам значение a.
  inline Vector_Nt& operator=(const T &a){
    for (int i=0; i<N; i++)
      v[i]=a;
    return *this;
  };

  //r сравнение. При отличии меньше чем на VEC_ZERO компоненты считаются одинаковыми
  inline bool operator==(const Vector_Nt &vect) const{
    for (int i=0; i<N ;i++)
      if(fabs(v[i]-vect.v[i])>VEC_ZERO)return false;
    return true;
  };

  inline bool operator!=(const Vector_Nt &vect) const{
    return (!(this->operator==(vect)));
  };

  inline Vector_Nt operator+(const Vector_Nt& vect) const{
    Vector_Nt result;
    for (int i=0; i<N ;i++)
      result.v[i]=v[i]+vect.v[i];
    return result;
  }

  inline Vector_Nt operator-(const Vector_Nt &vect) const {
    Vector_Nt result;
    for (int i=0; i<N ;i++)
      result.v[i]=v[i]-vect.v[i];
    return result;
  }
 
  //r Скалярное произведение векторов
  inline T operator*(const Vector_Nt& vect) const {
    T result=0;
    for (int i=0; i<N; i++)
      result+=v[i]*vect.v[i];
    return result;
  }

  //r Покомпонентное умножение на коэффициент
  inline Vector_Nt operator*(const T &coeff) const {
    Vector_Nt result;
    for (int i=0; i<N; i++)
      result[i]=coeff*v[i];
    return result;
  }

  //e vector multiplication (N=3 only)
  //r Векторное произведение
  inline Vector_Nt operator%(const Vector_Nt &r) const{ //reserved for N specializations
    if(N==3){
      return Vector_Nt(v[1]*r.v[2]-v[2]*r.v[1],v[2]*r.v[0]-v[0]*r.v[2],v[0]*r.v[1]-v[1]*r.v[0]);
    }
    return *this;
  }

  //r Умножение числа на вектор (переставлены местами множители).
 //  friend Vector_Nt operator*(T coeff,const Vector_Nt& vec);

  //r Покомпонентное деление на коэффициент
  inline Vector_Nt operator/(const T &coeff) const {
    Vector_Nt result;
    for (int i=0; i<N; i++)
      result[i]=v[i]/coeff;
    return result;
  }

  //r Умножение вектора на -1
  inline Vector_Nt operator-() const {
    Vector_Nt r;
    for (int i=0; i<N; i++)
      r.v[i]=-v[i];
    return r;
  }

  //r Сложение с присваиванием
  inline Vector_Nt& operator+=(const Vector_Nt &vect){
    for (int i=0; i<N; i++)
      v[i]+=vect.v[i];
    return *this;
  }

  //r Вычитание с присваиванием
  inline Vector_Nt& operator-=(const Vector_Nt &vect){
    for (int i=0; i<N; i++)
      v[i]-=vect.v[i];
    return *this;
  }

  //r Умножение на коэффициент с присваиванием
  inline Vector_Nt& operator*=(const T &coeff){
    for (int i=0; i<N; i++)
      v[i]*=coeff;
    return *this;
  }

  //r Деление на скаляр с присваиванием
  inline Vector_Nt& operator/=(const T &coeff){
    for (int i=0; i<N; i++)
      v[i]/=coeff;
    return *this;
  }

  //r Квадрат нормы вектора
  T norm2() const {
    T result=0;
    for (int i=0; i<N; i++)
      result+=v[i]*v[i];
    return result; 
  }

  //r Норма вектора
  T norm() const {
    return sqrt(norm2());
  }

  //r Возвращает норму и нормализует вектор (после этого его norm() вернет newnorm). 
  T normalize(T newnorm=1.){
    T norm=this->norm();
    if(norm>=VEC_ZERO){
      T c=newnorm/norm;
      for (int i=0; i<N; i++)
        v[i]*=c;
    }
    return norm;
  }

  //e Normalizes a vector that may have infinite components
  T inormalize(T newnorm=1.){
    T result=0;
    for (int i=0; i<N; i++){
      if(fabs(v[i])>=VEC_INFTY){
        if(result>=0)
          result=0.;
        result-=1;
        v[i]=v[i]>0 ? 1 : -1;
      }
      else if(result>=0) //still summing the normal components
        result+=v[i]*v[i];
      else
        v[i]=0.;
    }
    if(fabs(result)<VEC_ZERO)
      return 0.;
    newnorm/=sqrt(fabs(result));
    for (int i=0; i<N; i++)
      v[i]*=newnorm;
    return result<0 ? VEC_INFTY : result; 
  }


  //e nearest image distance within rectangular cell (FOR DISTANCE MEASUREMENTS)
  //e assumes that each coordinate absolute value is in the range [0,cell[i]) 
  //e returned vector is in the range [-cell[i]/2,cell[i]/2)
  //e flags indicate the periodicity in specific directions: 0x1 for X, 0x2 for Y, 0x4 for Z
  //r Ближайший образ в прямоугольной ячейке
  /*r Считаем, что все пространство разделено на ячейки - параллелепипеды с ребрами, параллельными
  осям координат и диагональю, заданной вектором rcell, причем начало координат является 
  центром одной из ячеек. Если *this находится центральной ячейке, возвращается копия *this.\n
  Иначе, если *this находится в кубе 3*3 ячейки с центром в начале координат, то создает образ 
  *this в центральной ячейке.\n
  Иначе, возвращает неопределенное значение. 
  */ 
  Vector_Nt rcell1(const Vector_Nt &cell,int flags=0xffff) const{
    Vector_Nt ret(*this);
    int i;
    for(i=0;i<N;i++){
      if(flags&(0x1<<i)){
        if(v[i]>cell[i]/2){
          ret[i]-=cell[i];
        }
        else if(v[i]<-cell[i]/2){
          ret[i]+=cell[i];
        }
      }
    }
    return ret;
  }

  //e @return a vector projection on a given axis
  Vector_Nt prj(int i) const {
    Vector_Nt res;
    res[i]=v[i];
    return res;
  }

  //e @return a vector projection on a given vector (k*v)*k
  Vector_Nt prj(const Vector_Nt &k) const {
    return k*(k*(*this));
  }

  //e @return a vector of length l parallel to this
  Vector_Nt unit(T l=1) const {
    Vector_Nt res(*this);
    res.normalize(l);
    return res;
  }

  //e reduction to elementary cell [0, cell[i]) (FOR REDUCTION TO ELEMENTARY CELL)
  //e flags indicate the periodicity in specific directions: 0x1 for X, 0x2 for Y, 0x4 for Z
  //r Почти то же, что и rcell1, но без ограничения на положение *this и с другой системой ячеек.
  /*r В начале координат находится не центр ячейки, а ее угол. Может работать медленнее из-за наличия 
  операции деления по модулю с плавающей точкой */ 
  Vector_Nt rcell(const Vector_Nt &cell, int flags=0xffff) const {
    Vector_Nt ret(*this);
    for (int i=0, flag=1; i<N; i++, flag<<=1) {
      if(flags&flag){
        ret.v[i]=fmod(v[i],cell[i]);
        if(ret.v[i]<0)ret.v[i]+=cell[i];
//        if(ret.v[i]<0 || ret.v[i]>cell[i])
  //        printf("!");
      }
    }
    return ret;
  }

  Vector_Nt rpcell(const Vector_Nt &p1, const Vector_Nt &cell, int flags=0xfff) const {
    Vector_Nt ret(*this);
    for (int i=0, flag=1; i<N; i++, flag<<=1) {
      if(flags&flag){
        if (ret.v[i]<p1[i] || ret.v[i]>p1[i]+cell[i]) {
          ret.v[i]=fmod(v[i]-p1[i],cell[i])+p1[i];
          if (ret.v[i]<p1[i]) ret.v[i]+=cell[i];
        }
      }
    }
    return ret;
  }
  

  //r Возвращает максимальную компоненту вектора  и ее индекс в ind 
  T maxcoord(int *ind=NULL) const {
    int im=0;
    T vv=v[0];
    for (int i=1; i<N; i++) {
      if(v[i]>vv){
        im=i;
        vv=v[i];
      }
    }
    if(ind)*ind=im;
    return vv;
  }


  //e returns the corrd having maximal absolute value
  T maxabscoord(int *ind=NULL) const {
    int im=0;
    T vv=fabs(v[0]);
    for (int i=1; i<N; i++) {
      if(fabs(v[i])>vv){
        im=i;
        vv=fabs(v[i]);
      }
    }
    if(ind)*ind=im;
    return v[im];
  }

  //e returns the corrd having minimal absolute value
  T minabscoord(int *ind=NULL) const {
    int im=0;
    T vv=fabs(v[0]);
    for (int i=1; i<N; i++) {
      if(fabs(v[i])<vv){
        im=i;
        vv=fabs(v[i]);
      }
    }
    if(ind)*ind=im;
    return v[im];
  }


  //r Возвращает минимальную компоненту вектора  и ее индекс в ind 
  T mincoord(int *ind=NULL) const {
    int im=0;
    T vv=v[0];
    for (int i=1; i<N; i++) {
      if(v[i]<vv){
        im=i;
        vv=v[i];
      }
    }
    if(ind)*ind=im;
    return vv;
  }

  //r Выводит вектор в поток вывода по умолчанию, в формате (x,y,z)\\n
  void print() const{
    cout<< "(";
    for(int i=0;i<N;i++){
      cout<< v[i];
      if(i!=N-1)
        cout<< ", ";
    }
    cout<< ")\n";
  }

  //e returns true if the vector has infinite components
  bool infinite() const {
    for(int i=0;i<N;i++){
      if(fabs(v[i])>=VEC_INFTY)
        return true;
    }
    return false;
  }
};

template<class T , int N>
Vector_Nt<T, N> operator*(const T &coeff,const Vector_Nt<T, N> &vec){
  return vec*coeff;
}

template<class T , int N>
Vector_Nt<T, N> operator*(int coeff,const Vector_Nt<T, N> &vec){
  return vec*coeff;
}

//template <> Vector_Nt<double, 3> operator*(const double &coeff,const Vector_Nt<double,3> &vec);

// old Vector_3 compatibility typedefs and functions
typedef Vector_Nt<vec_type, 3> Vector_3;
typedef Vector_3 *Vector_3P;
typedef Vector_Nt<vec_type, 2> Vector_2;
typedef Vector_2 *Vector_2P;

template <int N> 
class  Vector_N: public Vector_Nt<vec_type, N>{
};

//Vector_3 operator*(int coeff,const Vector_3 &vec){
//  return vec*((vec_type)coeff);
//}


//e finds the maximum distance between vector pairs
//r Находит максимальное расстояние между векторами va1[i], va2[i], i=1..n
/*r @param va1 - массив Vector_3[n]
    @param n - длина массивов va1 и va2
*/
vec_type dist_max(Vector_3 *va1,Vector_3 *va2,int n);

//e finds average distance between vector pairs
//r Находит среднее расстояние между векторами va1[i], va2[i], i=1..n
vec_type dist_av(Vector_3 *va1,Vector_3 *va2,int n);

//e finds the average difference norm between two vector sets of the same length
/*e optionally gives the indexes for maximal and minimal difference
 va2 can be NULL, then the norm of va1 is used */

//r Находит среднее расстояние между va1[i] и va2[i], а также, по желанию, индексы, на которых достигается min и max расстояние
vec_type diff_av(Vector_3 *va1,Vector_3 *va2,int n, int *minind=0, int *maxind=0);

//e finds suitable perpendicular to a vector
//r Находит перпендикуляр к вектору vAB
Vector_3 FindPerp(const Vector_3 &vAB);



//e Returns the average (center) vector of the vector array
//e and cooordinates of a minimal box in which it is contained
Vector_3 GetScope(const Vector_3 *varr,long n,Vector_3* box_min,Vector_3* box_max);

//e the same with long index array
Vector_3 GetIScope(const Vector_3 *varr,long *indarr,long n,Vector_3* box_min,Vector_3* box_max);

//e the same with int index array
Vector_3 GetIScopei(const Vector_3 *varr,int *indarr,int n,Vector_3* box_min,Vector_3* box_max);

// neue Funktionen

//e clears vector array with optional integer index
//r Очистка массива векторов, с поддержкой индексирования 
/*r 
В данном Vector_3 vec[] обнуляет n координат. Если ind==NULL, то 
очищает первые n элементов. Если ind!=NULL, то для i=0..n-1
очищает vec[ind[i]]
См. @ref indexed_calculations.
*/ 
void clear_vecarri(int n,Vector_3 *vec, int *ind=0);

//e reflects the vector ini+dir*t+0.5*force*t^2 to be inside a box limited by 0 and box sizes
//e changes dir according to the final state
//e fills crossed dir with bit flags corresponding directions along which the walls were crossed
Vector_3 Reflect(Vector_3& ini, double t,Vector_3 &dir, double *box, int flag=0x7, const Vector_3 &force=Vector_3()); 


inline vec_type vec_area(const Vector_2 &vect1, const Vector_2 &vect2) {
  return fabs(vect1[0]*vect2[1]-vect1[1]*vect2[0])/2;
};

inline vec_type vec_area(const Vector_3 &vect1, const Vector_3 &vect2) {
  return (vect1%vect2).norm()/2;
};

// remake for vec_type!

template <class num_t, class denum_t, class val_t>
denum_t acccomp(const num_t x, const denum_t y, const val_t val) {
  num_t eps_num=numeric_limits<num_t>::epsilon();
  denum_t eps_denum=numeric_limits<denum_t>::epsilon();
  return (eps_num>=eps_denum) ? acccomp(x, (num_t)y, (num_t)val) : acccomp((denum_t)x, y, (denum_t)val);
}

template <class num_t, class denum_t>
denum_t acccomp(const num_t x, const denum_t y) {
  return acccomp(x, y, num_t(1));
}

template<class T>
bool acccomp(const T x, const T y, const T val) {
  T eps=numeric_limits<T>::epsilon();
  int iexp;
  frexp(val,&iexp);
  T mult=ldexp(.5, iexp+1);
  T err=T(256)*mult*eps;
  return fabs(x-y)<=err;
}

template<class T>
bool acccomp(const T x, const T y) {
  return acccomp(x, y, T(1));
}

template <class num_t, class denum_t>
denum_t accdiv(const num_t x, const denum_t y) {
  num_t eps_num=numeric_limits<num_t>::epsilon();
  denum_t eps_denum=numeric_limits<denum_t>::epsilon();
  return (eps_num>=eps_denum) ? accdiv1(x, (num_t)y) : accdiv1((denum_t)x, y);
}


template <class T>
T accdiv1(T x, T y) {
  T result;
  T eps=numeric_limits<T>::epsilon();

  T fr=x/y;
  T ifr=floor(fr+T(.5));
//  T err=64*eps*(ifr!=0 ? fabs(ifr) : 1);
  int mult;
  if (ifr<=512)
    mult=512;
  else {
    int iexp;
    frexp(ifr,&iexp);
    mult=int(ldexp(.5, iexp+1));
  }
  T err=mult*eps;
  if (fabs(fr-ifr)<=err)
    result=ifr;
  else
    result=fr;
  return result;
}




template <class num_t, class denum_t, class P>
denum_t accdiv_rmd(const num_t x, const denum_t y, P *remainder) {
  num_t eps_num=numeric_limits<num_t>::epsilon();
  denum_t eps_denum=numeric_limits<denum_t>::epsilon();
  return (eps_num>=eps_denum) ? accdiv_rmd(x, (num_t)y, remainder) : accdiv_rmd((denum_t)x, y, remainder);
}

template <class T, class P>
T accdiv_rmd(const T x, const T y, P *remainder) {
  T result=accdiv(x, y);
  T iresult=floor(result);
  if (result-iresult>0)
    *remainder=x-iresult*y;
  else
    *remainder=0;
  return result;
}

//e returns random unit vector uniformely distributed in space (?? check this)
inline Vector_3 randdir(){
  vec_type xi1=2*((vec_type)rand())/RAND_MAX-1.;
  vec_type xi2=((vec_type)rand())/RAND_MAX;
  vec_type r1=sqrt(1.-xi1*xi1);
  return  Vector_3(r1*cos(2*M_PI*xi2),r1*sin(2*M_PI*xi2),xi1);
}

///\en Calculates extent of the vector container.
///    \return the center of the vector set, optionally
///    (if arguments are not NULL) fills the bounding box in \a box_min, \a box_max.
template<class vec_inp_it>
Vector_3 get_extent(vec_inp_it beg,vec_inp_it end, Vector_3* box_min=NULL,Vector_3* box_max=NULL){
  if(beg==end)
    return Vector_3();
  Vector_3 center(*beg++);
  Vector_3 cube1(center), cube2(center);
  size_t n=1;
  for(;beg!=end;++beg){
    Vector_3 vec=*beg;
    center+=vec;
    for(size_t j=0;j<3;j++){
      if(cube1[j]>vec[j])
        cube1[j]=vec[j];
     if(cube2[j]<vec[j])
        cube2[j]=vec[j];
    }
    n++;
  }
  if(box_min)
    *box_min=cube1;
  if(box_max)
    *box_max=cube2;
  return center/n;
}

///\en Returns \a true if absolute value of \a a is at least \a factor times less than
///    the absolute value of \a b.
template<class T>
bool much_less(const T& a, const T& b, const T &factor){
  if(fabs(a)<fabs(b))
    return fabs(a*factor/b)<1 ? true : false;
  else
    return false;
}


# endif // __VECTOR_3_H

