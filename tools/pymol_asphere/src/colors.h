/***************************************************************************
                                  colors.h
                              W. Michael Brown
                             -------------------

  Storage and manipulation of colors.                             

    begin                : Wed Jun 11 2003
    copyright            : (C) 2003 by mbrown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#ifndef COLORS_H
#define COLORS_H

#include "cartesian.h"
#include <map>
#include <string>
#include <iostream>
using namespace std;

/// Class for dealing with RGB color space
class Colors {
public: 
	Colors();
	~Colors();

  /// Return number of colors in the map
  unsigned size();
  
  /// Return the color in the gradient between start and end 
  /** \param start Starting value
    * \param end Ending Value
    * \param value Color in gradient determined for value
    * relative to start and end
    * \param startcolor Starting color for the gradient
    * \param endcolor Ending color for the gradient */
  colorPt gradient(double start, double end, double value,
                   const colorPt &startcolor, const colorPt &endcolor);
  /// Three color gradient
  /** \sa gradient() */
  colorPt gradient(double start, double mid, double end, double value,
                   const colorPt &startcolor, const colorPt &midcolor,
                   const colorPt &endcolor);

  /// Output the colors to standard out
  void outputcolors();
  /// Return a string with the list of available colors
  string colorlist();
  
  /// Return an RGB color based on name
  colorPt operator[](const string &name);
  
private:
  map<string, colorPt> colors;
};

#endif

