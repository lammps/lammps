#ifndef STRINGMANIP_H
#define STRINGMANIP_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>

using std::cout;
using std::transform;
using std::fstream;
using std::vector;
using std::string;

//* A collection of string convesion/manipulation routines 
namespace ATC_STRING
{
  //* converts anything that has iostream::<< defined to a string 
  template<typename T>
  inline string tostring(const T &v, int precision=0)
  {
    std::ostringstream out;
    if (precision) out << std::setprecision(precision);
    out << v;
    return out.str();
  }

  //* convert a string to anything that has iostream::>> defined
  template<typename T>
  inline T str2T(const string &s, T v)
  {
    std::istringstream in(s);
    if (!(in >> v)) cout<<"Error: str2T invalid string conversion\n";\
    return v;
  }

  //* convert a string to a double
  inline double str2dbl(const string &s) { return str2T(s, double(0.0)); }
  
  //* convert a string to an int
  inline int str2int(const string &s)    { return str2T(s, int(0)); }

  //* replaces all characters in a set of characters with a character
  //* @param str input string
  //* @param *s pointer to array of characters that will be replaced
  //* @param r character to replace characters in s[] with
  static void replace_all_of(string &str, const char *s, const char r)
  {
     int found;
     found = str.find_first_of(s);
     while (found != string::npos)
     {
        str[found] = r;
        found = str.find_first_of(s, found+1);
     }
  }

  //* converts the string to lowercase
  inline string& to_lower(string &s)
  {
    transform(s.begin(),s.end(),s.begin(),static_cast<int(*)(int)>(tolower));
    return s;
  }

  //* converts the string to uppercase
  inline string& to_upper(string &s)
  {
    transform(s.begin(),s.end(),s.begin(),static_cast<int(*)(int)>(toupper));
    return s;
  }
  
  //* removes any whitespace from the beginning or end of string
  static string& trim(string &s)
  {
    if (s.empty()) return s;
    size_t found = s.find_last_not_of(" \t\n");
    if (found < s.size()-1) s.erase(found+1);
    found = s.find_first_not_of(" \t\n");
    if (found != string::npos) s = s.substr(found);
    return s;
  }

  //* splits delimited string into a vector of strings
  static void split_up_string(const string &s, vector<string> &ss, char del=' ')
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
  
  //* scans a string for a list of commands delimited by ', \t' with # comments
  //* @param line The input line to be parsed
  //* @cs A vector of strings parsed from the input line
  static void get_command_strings(string line, vector<string> &cs)
  {
    if (line.find('#') != string::npos)  line.erase(line.find("#"));
    replace_all_of(line, "\t,", ' ');
    to_lower(trim(line));
    split_up_string(line, cs);
  }
 
  //* reads a single line from a file and splits it into a vector of strings
  //* returns the number of strings in the vector
  inline int get_command_line(fstream &fid, vector<string> &cs)
  {
     string line;
     getline(fid, line);
     get_command_strings(line, cs);
     return cs.size();
  }
}

#endif
