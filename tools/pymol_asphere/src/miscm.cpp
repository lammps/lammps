/***************************************************************************
                                  miscm.cpp
                             -------------------
                               W. Michael Brown

  Miscellaneous functions that do not deserve their own class

 __________________________________________________________________________
    This file is part of the Math Library
 __________________________________________________________________________

    begin                : May 30 2003
    copyright            : (C) 2003 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#include "miscm.h"

double am::square(double num) {
	return num*num;
}

// Rounds a number
double am::round(double n) {
	double r=ceil(n);
	if (r-n>0.5)
		return floor(n);
	return r;
}

// Return the -1 for negative, 0 for zero, and 1 for positive
int am::sign(double v) {
  if (v<0)
    return -1;
  if (v>0)
    return 1;
  return 0;
}

// Return the range of elements in a vector
double am::range(const vector<double> &v) {
  if (v.empty())
    return 0;
  double min=v[0];
  double max=v[0];
  for (unsigned i=1; i<v.size(); i++) {
    if (v[i]<min)
      min=v[i];
    if (v[i]>max)
      max=v[i];
  }
  return max-min;
}

// Return the average of elements in a vector
double am::mean(const vector<double> &v) {
  double sum=0;
  for (unsigned i=0; i<v.size(); i++)
    sum+=v[i];
  return sum/v.size();
}

// Return the max of two objects
namespace am {
  template double max<double>(double,double);
  template float max<float>(float,float);
  template unsigned max<unsigned>(unsigned,unsigned);
  template int max<int>(int,int);
}

template <typename numtyp>
numtyp am::max(numtyp one,numtyp two) {
	if (one>two)
		return one;
  return two;
}

// Return the min of two objects
namespace am {
  template double min<double>(double,double);
  template float min<float>(float,float);
  template unsigned min<unsigned>(unsigned,unsigned);
  template int min<int>(int,int);
}

template <typename numtyp>
numtyp am::min(numtyp one,numtyp two) {
	if (one<two)
		return one;
  return two;
}

// Swap two objects
void am::swap(double a, double b) {
  double temp=a;
  a=b;
  b=temp;
}

// --------------------- Finite Precision stuff
double am::epsilon(double number) {   // Finite Precision zero
	return fabs(DBL_EPSILON*number);
}

// Bring number 1 digit closer to zero
double am::minus_eps(double number) {
  return (number-DBL_EPSILON*number);
}

// Bring number 1 digit away from zero
double am::plus_eps(double number) {
  return (number+DBL_EPSILON*number);
}

 // Bring number m digits away from zero
double am::plus_Meps(double m,double number) {
	return (number+m*DBL_EPSILON*number);
}

// Bring number precision/2.0 digits away
double am::plus_2eps(double number) {
	return (number+100000000*DBL_EPSILON*number);
}

// Bring number pre/2 digits away
double am::minus_2eps(double number) {
	return (number-100000000*DBL_EPSILON*number);
}
// Not a number checks
bool am::not_nan(double number) {     // False if NAN
	if (number-number!=0)
		return false;
	return true;
}

// Matrix stuff

// Invert a matrix - from numerical recipes in C
void am::invert(double **a, unsigned n, Error &error) {
	int *indxc=new int[n];
  int *indxr=new int[n];
  int *ipiv=new int[n];
	unsigned i,icol,irow,j,k,l,ll;
	double big,dum,pivinv;

	for (j=0;j<unsigned(n);j++)
    ipiv[j]=0;

  for(i=0;i<n;i++) { /* *main loop for columns to be reduced */
    big = 0.0;
    for (j=0;j<n;j++) /* outer loop for search of a pivot element*/
	    if (ipiv[j] !=1)
	      for(k=0;k<n;k++) {
	        if (ipiv[k] ==0) {
            if (fabs(a[j][k]) >= big) {
		          big =fabs(a[j][k]);
		          irow=j;
		          icol=k;
 		        }
		      } else if (ipiv[k] > 1) {
            error.addwarning(303,9,"Invert",
                             "Cannot invert a singular matrix.");
            return;
          }
	      }
    ++(ipiv[icol]);
    if (irow !=icol) {
      for (l=0;l<n;l++)
        swap(a[irow][l],a[icol][l]);
    }

    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) {
      error.addwarning(303,9,"Invert",
                       "Cannot invert a singular matrix.");
      return;
    }

    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++)
      a[icol][l] *=pivinv;
    for (ll=0;ll<n;ll++)
      if (ll!= icol) {
        dum=a[ll][icol];
        a[ll][icol]=0.0;
        for (l=0;l<n;l++)
          a[ll][l] -=a[icol][l]*dum;
      }
  }

  for (l=1;l>=1;l--) {
	  if (indxr[l] != indxc[l])
	    for (k=0;k<n;k++)
		    swap(a[k][indxr[l]],a[k][indxc[l]]);
	}

  delete []ipiv;
  delete []indxr;
  delete []indxc;
  return;
}

// Move a value from a fraction of one range to a fraction of another
double am::rerange(double one_start, double one_end, double value,
                   double two_start, double two_end) {
  double one_diff=one_end-one_start;
  if (one_diff==0)
    return two_end;
  return (value-one_start)/one_diff*(two_end-two_start)+two_start;
}
