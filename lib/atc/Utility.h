#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <ctime>
#include <limits>
#undef near

#include <vector>
using std::vector;
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "math.h"

#include "ATC_Error.h"

using std::transform;
using std::fstream;
using std::vector;
using std::string;
using std::strcmp;


namespace ATC_Utility
{
  /** constants */
  static const double Pi = 4.0*atan(1.0);
  static const double Big = 1.e20;
  const static double parsetol = 1.0e-10;  
  const static double parsebig = 1.0e10;

  /** scalar triple product */
  inline double det3(double * v1, double * v2, double * v3) {
    return
     -v1[2]*v2[1]*v3[0] + v1[1]*v2[2]*v3[0] + v1[2]*v2[0]*v3[1]
     -v1[0]*v2[2]*v3[1] - v1[1]*v2[0]*v3[2] + v1[0]*v2[1]*v3[2];
  }
  inline void plane_coords(int i, int & i1, int & i2) {
    if      (i==0) { i1 = 1; i2 = 2; }
    else if (i==1) { i1 = 2; i2 = 0; }
    else if (i==2) { i1 = 0; i2 = 1; }
  }
  /** 2d cross product */
  inline double cross2(double * v1, double * v2) {
    return v1[0]*v2[1] - v1[1]*v2[0];
  }
  /** 3d cross product */
  inline void cross3(double * v1, double * v2, double * v) {
    v[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v[1] = v1[2]*v2[0] - v1[0]*v2[3];
    v[2] = v1[0]*v2[1] - v1[1]*v2[0];
  }
  inline double norm3(double * v) {return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }

  inline int sgn(double x) { return (int) ((x>0) - (x<0)); }
  inline int rnd(double x) { return (int) (x+sgn(x)*0.5); }

  /** Compares doubles acceptably */
  inline bool dbl_geq(double dblL, double dblR, int tolMult=2)
  {
    return ( (dblL >   dblR) ||
            ((dblL <= (dblR + \
                       std::numeric_limits<double>::epsilon() * \
                       tolMult * std::max(abs(dblL), abs(dblR)))) && 
             (dblR <= (dblL + \
                       std::numeric_limits<double>::epsilon() * \
                       tolMult * std::max(abs(dblL), abs(dblR))))));
  }

    
  inline double parse_min(const char * arg) {
    if (strcmp(arg,"INF") == 0)  return -parsebig;
    else                         return atof(arg);  }   
  inline double parse_max(const char * arg) {    
     if (strcmp(arg,"INF") == 0)  return parsebig;    
     else                         return atof(arg);
  }   
  inline double parse_minmax(const char * arg) {
    if      (strcmp(arg,"-INF") == 0)  return -parsebig;
    else if (strcmp(arg,"INF") == 0)   return  parsebig;
    else                         return atof(arg);  }   
  inline void split_values(double & min, double & max) {
    double eps = std::max(parsetol,parsetol*fabs(min));    
    min -= eps;
    max += eps;
  }   
    

  /** Returns true if the value v is between min &  max */
  template<typename T>
  inline bool in_range(const T &v, const T &min, const T &max)
  {
    return v >= min && v <= max;
  }
  /** Returns true if the value v is between min & max within tolerance TOL */
  template<typename T>
  inline bool in_range(const T &v, const T &min, const T &max, const T &tol)
  {
    return in_range(v, min-tol, max+tol);
  }
  /** Returns the value with the larger absolute value */
  template<typename T>
  inline T max_abs(const T &a, const T &b)
  {
    return (a*a > b*b) ? a : b;
  }
  /** Returns the value with the smaller absolute value */
  template<typename T>
  inline T min_abs(const T &a, const T &b)
  {
    return (a*a < b*b) ? a : b;
  }
  /** A simple Matlab like timing function */
  inline double tictoc()
  {
    double t_new = clock() / (double) CLOCKS_PER_SEC;
    static double t = t_new;
    double dt = t_new - t;
    t = t_new;
    return dt; 
  }
  /** A simple timer */
  inline double timer()
  {
    double t_new = clock() / (double) CLOCKS_PER_SEC;
    static double t = t_new; // done once at first time the function is evoked
    double dt = t_new - t;
    return dt; 
  }
  /** Binary search between low & high for value (assumes array is sorted) */
  
  template<typename T>
  inline int search_sorted(const T* A, T value, int low, int high)
  {
    int mid;
    while (low < high)
    {
     mid = (low + high) >> 1;                 // set mid to mean of high and low, ">>" is a bit shift which divides by 2, rounding down
     if (A[mid] > value)       high = mid;  // mid was too high, reduce high
     else if (A[mid] < value)  low = mid+1;   // mid was too low, increase low
     else return mid;
    }
    return -1; // value not found in array, return -1
  }
  /** Flat search between low & high for value (assumes array is NOT sorted) */
  template<typename T>
  inline int search_unsorted(const T* A, T value, int low, int high)
  {
    for (int i=low; i<high; i++) if (A[i] == value) return i;
    return -1;
  }
  /** Regular swap */
  template<typename T>
  inline void swap(T &a, T &b)
  {
    T temp = a;
    a = b;
    b = temp;
  }


  //===================================================================
  /** A collection of string convesion/manipulation routines  */
  //===================================================================

  /** converts anything that has iostream::<< defined to a string  */
  template<typename T>
  inline string to_string(const T &v, int precision=0)
  {
    std::ostringstream out;
    if (precision) out << std::setprecision(precision);
    out << v;
    return out.str();
  }

  /** conversion to string */
  inline string true_false(const double &v)
  {
    if (v) return "TRUE";
    else   return "FALSE";
  }

  /** test conversion to double */
  inline bool is_numeric(const string &s)
  {
    double v;
    std::istringstream in(s);
    return (in >> v);
  }

  /** convert a string to anything that has iostream::>> defined */
  template<typename T>
  inline T str2T(const string &s, T v) 
  {
    std::istringstream in(s);
    if (!(in >> v)) throw ATC::ATC_Error("str2T invalid string conversion");
    return v;
  }

  /** convert a string to a double */
  inline double str2dbl(const string &s) { return str2T(s, double(0.0)); }

  /** tests if conversion to double is possible */
  inline bool is_dbl(const string &s) 
  {
    char *endptr;
    strtod(s.c_str(), &endptr);
    if(endptr != NULL && *endptr == '\0') return true;
    return false;
  }
  
  /** convert a string to an int */
  inline int str2int(const string &s)    { return str2T(s, int(0)); }

  //* replaces all characters in a set of characters with a character
  //* @param str input string
  //* @param *s pointer to array of characters that will be replaced
  //* @param r character to replace characters in s[] with
  static void replace_all_of(string &str, const char *s, const char r)
  {
     size_t found;
     found = str.find_first_of(s);
     while (found != string::npos)
     {
        str[found] = r;
        found = str.find_first_of(s, found+1);
     }
  }

  /** converts the string to lowercase */
  inline string& to_lower(string &s)
  {
    transform(s.begin(),s.end(),s.begin(),static_cast<int(*)(int)>(tolower));
    return s;
  }

  /** converts the string to uppercase */
  inline string& to_upper(string &s)
  {
    transform(s.begin(),s.end(),s.begin(),static_cast<int(*)(int)>(toupper));
    return s;
  }
  
  /** removes any whitespace from the beginning or end of string */
  static string& trim(string &s)
  {
    if (s.empty()) return s;
    size_t found = s.find_last_not_of(" \t\n");
    if (found < s.size()-1) s.erase(found+1);
    found = s.find_first_not_of(" \t\n");
    if (found != string::npos) s = s.substr(found);
    return s;
  }

  /** splits delimited string into a vector of strings */
  static void split(const string &s, vector<string> &ss, char del=' ')
  {
    ss.clear();
    size_t begin=0, end=0;
    while (end != string::npos)
    {
      begin = s.find_first_not_of(del, end);     // find beginning of fragment
      end = s.find_first_of(del,begin);
      if (begin != string::npos)                 // safe if end is npos-1 
        ss.push_back(s.substr(begin,end-begin)); 
    }
  } 

  static string cap(string & s) {
    s[0] = toupper(s[0]);
    return s;
  }

  /* turns Aa_Bb_Cc into aa_bb_cc */
  static string to_cap(const string &s) {
    vector<string> words;
    string delimiter = "_";
    split(s,words,(delimiter.c_str())[0]);
    string name = "";
    for (unsigned int i = 0; i < words.size(); i++) {
      name = name + cap(words[i]);
    }
    return name;
  }
  
  //* scans a string for a list of commands delimited by ', \t' with # comments
  //* @param line The input line to be parsed
  //* @cs A vector of strings parsed from the input line
  static void command_strings(string line, vector<string> &cs)
  {
    if (line.find('#') != string::npos)  line.erase(line.find("#"));
    replace_all_of(line, "\t,", ' ');
    to_lower(trim(line));
    split(line, cs);
  }
 
  //* reads a single line from a file and splits it into a vector of strings
  //* returns the number of strings in the vector
  inline int command_line(fstream &fid, vector<string> &cs)
  {
     string line;
     getline(fid, line);
     command_strings(line, cs);
     return cs.size();
  }
}

#endif
