/*
 Atz_VTK_Helper.h
 
 Paul J. Atzberger
 http://atzberger.org/
 
 
 Collection of routines used for help
 in parsing the VTK data sets.
 
*/

#ifndef ATZ_VTK_HELPER_H_
#define ATZ_VTK_HELPER_H_


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "SELM_Package.h"

using namespace std; /* ensures standard template library in namespace */
//using namespace Atz_VTK_Helper; /* ensures standard template library in namespace */
using namespace LAMMPS_NS;

class Atz_VTK_Helper  {

public:

  static const char *error_str_code;

  static int writeVecFieldVTK_XML_File(const char   *filename,
                                       int           num_dim,
                                       int          *numMeshPtsPerDir,
                                       double       *meshCenterX0,
                                       double        meshDeltaX,
                                       const char   *vec_name,
                                       double      **vec_array);

  static int writeScalarFieldVTK_XML_File(const char   *filename,
                                          int           num_dim,
                                          int          *numMeshPtsPerDir,
                                          double       *meshCenterX0,
                                          double        meshDeltaX,
                                          const char   *scalar_name,
                                          double       *scalar_array);
protected:

  public:
  Atz_VTK_Helper();
  virtual ~Atz_VTK_Helper();

  protected:

  public:


};

#endif /* ATZ_VTK_HELPER_H_ */

