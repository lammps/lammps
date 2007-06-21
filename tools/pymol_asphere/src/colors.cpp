/***************************************************************************
                                 colors.cpp
                              W. Michael Brown
                             -------------------

    begin                : Wed Jun 11 2003
    copyright            : (C) 2003 by mbrown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#include "colors.h"

Colors::Colors() {
  colors.insert(pair<string,colorPt>("white",colorPt(1,1,1)));
  colors.insert(pair<string,colorPt>("red",colorPt(1,0,0)));
  colors.insert(pair<string,colorPt>("green",colorPt(0,1,0)));
  colors.insert(pair<string,colorPt>("blue",colorPt(0,0,1)));
  colors.insert(pair<string,colorPt>("cmyk_blue",colorPt(0.074,0.168,0.614)));
  colors.insert(pair<string,colorPt>("yellow",colorPt(1,1,0)));
  colors.insert(pair<string,colorPt>("purple",colorPt(0.697,0.2,0.697)));
  colors.insert(pair<string,colorPt>("hotpink",colorPt(1,0,0.5)));
  colors.insert(pair<string,colorPt>("brown",colorPt(0.551,0.272,0.15)));
  colors.insert(pair<string,colorPt>("orange",colorPt(1,0.5,0)));
  colors.insert(pair<string,colorPt>("black",colorPt(0,0,0)));
  colors.insert(pair<string,colorPt>("magenta",colorPt(1,0,1)));
  colors.insert(pair<string,colorPt>("slate",colorPt(0.5,0.5,1)));
  colors.insert(pair<string,colorPt>("marine",colorPt(0,0.5,1)));
  colors.insert(pair<string,colorPt>("cmyk_marine",colorPt(0.273,0.469,0.758)));
  colors.insert(pair<string,colorPt>("teal",colorPt(0.2,0,0.598)));
  colors.insert(pair<string,colorPt>("forest",colorPt(0.098,0.5,0.098)));
  colors.insert(pair<string,colorPt>("deep",colorPt(0.098,0.098,0.5)));
  colors.insert(pair<string,colorPt>("grey",colorPt(0.5,0.5,0.5)));
  colors.insert(pair<string,colorPt>("wheat",colorPt(0.988,0.819,0.646)));
}

Colors::~Colors(){
}

// Number of colors in the map
unsigned Colors::size() {
  return colors.size();
}  

// Return the color in the gradient between start and end according to
//   value
colorPt Colors::gradient(double start, double end, double value,
                         const colorPt &startcolor, const colorPt &endcolor) {
  double percent;

  if (value<start)
    return startcolor;
  else if (value>end)
    return endcolor;
  if (start==end)
    return startcolor;
    
  percent=((value-start)/(end-start));
    return ((endcolor-startcolor)*percent)+startcolor;
}
  
// Three color gradient
colorPt Colors::gradient(double start, double mid, double end, double value,
                         const colorPt &startcolor, const colorPt &midcolor,
                         const colorPt &endcolor) {
  double percent;

  if (value==mid)
    return midcolor;
  else if (value<start)
    return startcolor;
  else if (value>end)
    return endcolor;
  else if (value<mid) {
    percent=((value-start)/(mid-start));
    return ((midcolor-startcolor)*percent)+startcolor;
  } else {
    if (mid==end)
      return midcolor;
    percent=((value-mid)/(end-mid));
    return ((endcolor-midcolor)*percent)+midcolor;
  }
}

// Output the colors to standard out
void Colors::outputcolors() {
    map<string, colorPt>::iterator m;

    for (m=colors.begin(); m!=colors.end(); m++)
      cout << "    " << m->first << endl;
}
 
// Return a string with the list of available colors
string Colors::colorlist() {
  string colorl;
  map<string, colorPt>::iterator m;

  for (m=colors.begin(); m!=colors.end(); m++)
    colorl+="\n\t"+m->first;
  colorl+="\n";
  return colorl;
}

colorPt Colors::operator[](const string &name) {
    map<string, colorPt>::iterator m;
    m=colors.find(name);
    if (m==colors.end())
      return colorPt(1,1,1); // Return white if not found
    else
      return m->second;
}

