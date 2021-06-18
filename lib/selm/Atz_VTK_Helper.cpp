/*

  Atz_VTK_Helper.cpp 
 
  Paul J. Atzberger
  http://atzberger.org/
 
*/

#include "Atz_VTK_Helper.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Atz_XML_Package.h"

const char *Atz_VTK_Helper::error_str_code = "Atz_VTK_Helper.cpp";

Atz_VTK_Helper::Atz_VTK_Helper() {
  // TODO Auto-generated constructor stub
}


Atz_VTK_Helper::~Atz_VTK_Helper() {
  // TODO Auto-generated destructor stub
}


int Atz_VTK_Helper::writeVecFieldVTK_XML_File(const char   *filename,
                                              int           num_dim,
                                              int          *numMeshPtsPerDir,
                                              double       *meshCenterX0,
                                              double        meshDeltaX,
                                              const char   *vec_name,
                                              double      **vec_array) {

  const char *error_str_func = "writeVecFieldVTK_XML_File()";

  stringstream extrasStr;
  stringstream nameTag;

  ofstream fid;

  // do some precalculations
  int N  = numMeshPtsPerDir[0]*numMeshPtsPerDir[1]*numMeshPtsPerDir[2];
  int Nw = numMeshPtsPerDir[0]*numMeshPtsPerDir[1]*(numMeshPtsPerDir[2]-1);

  int numTotalMeshPts = N;

  int J[3];
  int index;

  double deltaX;
  double X[3];

  // open the file
  fid.open(filename);

  if (!fid.is_open()) {  // if file is not open
    stringstream message;
    message << "Could not open file to write error ocurred." << endl;
    message << "  filename = " << filename << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  Atz_XML_Package::writeXMLHeader(fid);

  nameTag.str("VTKFile");
  extrasStr.str("");
  extrasStr << "type=\""       << "RectilinearGrid" << "\" ";
  extrasStr << "version=\""    << "0.1" << "\" ";
  extrasStr << "byte_order=\"" << "LittleEndian" << "\"";
  Atz_XML_Package::writeTagStart(fid,
                                 nameTag.str().c_str(),
                                 extrasStr.str().c_str());

  nameTag.str("RectilinearGrid");
  extrasStr.str("");
  extrasStr << "WholeExtent=\"" << "0 " << (numMeshPtsPerDir[0] - 1) << " "
                                << "0 " << (numMeshPtsPerDir[1] - 1) << " "
                                << "0 " << (numMeshPtsPerDir[2] - 1) << "\" ";
  extrasStr << "GhostLevel=\""  << "0" << "\" ";
  Atz_XML_Package::writeTagStart(fid,
                                 nameTag.str().c_str(),
                                 extrasStr.str().c_str());

  nameTag.str("Piece");
  extrasStr.str("");
  extrasStr << "Extent=\"" << "0 " << (numMeshPtsPerDir[0] - 1) << " "
                           << "0 " << (numMeshPtsPerDir[1] - 1) << " "
                           << "0 " << (numMeshPtsPerDir[2] - 1) << "\" ";
  Atz_XML_Package::writeTagStart(fid,
                                 nameTag.str().c_str(),
                                 extrasStr.str().c_str());

  nameTag.str("Coordinates");
  extrasStr.str("");
  //extrasStr << "Scalars=\"" << "(Name of default scalar field here)" << "\" ";
  extrasStr << "Vectors=\"" << vec_name << "\"";
  Atz_XML_Package::writeTagStart(fid,
                                 nameTag.str().c_str(),
                                 extrasStr.str().c_str());

  // points along the x-axis, y-axis, z-axis
  for (int d = 0; d < num_dim; d++) {
    double L       = numMeshPtsPerDir[d]*meshDeltaX;
    double baseX_d = (meshCenterX0[d] - 0.5*L) + 0.5*meshDeltaX;
    //double baseX_d = (meshCenterX0[d] - 0.5*L);
    double X_d;
    double RangeMin = baseX_d;
    double RangeMax = baseX_d + (numMeshPtsPerDir[d] - 1)*meshDeltaX;

    nameTag.str("DataArray");
    extrasStr.str("");
    extrasStr << "type=\"" << "Float32" << "\" ";
    switch (d) {
    case 0:
      extrasStr << "Name=\"" << "x" << "\" ";
      break;
    case 1:
      extrasStr << "Name=\"" << "y" << "\" ";
      break;
    case 2:
      extrasStr << "Name=\"" << "z" << "\" ";
      break;
    default:
      extrasStr << "Name=\"" << "unknown-axis" << "\" ";
      break;
    } // end switch
    extrasStr << "format=\"" << "ascii" << "\" ";
    extrasStr << "RangeMin=\"" << RangeMin << "\" ";
    extrasStr << "RangeMax=\"" << RangeMax << "\"";
    Atz_XML_Package::writeTagStart(fid, nameTag.str().c_str(),
                                   extrasStr.str().c_str());

    // write the d-axis coordinates
    for (int m = 0; m < numMeshPtsPerDir[d]; m++) {
      X_d = baseX_d + m*meshDeltaX;
      fid << X_d << " ";
    }
    fid << endl;

    nameTag.str("DataArray");
    Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  } // end of the d-loop over different directions

  nameTag.str("Coordinates");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  nameTag.str("PointData");
  Atz_XML_Package::writeTagStart(fid,nameTag.str().c_str());
  // consider now each output array and if the flag is set

  nameTag.str("DataArray");
  extrasStr.str("");
  extrasStr << "type=\"" << "Float32" << "\" ";
  extrasStr << "Name=\"" << vec_name << "\" ";
  extrasStr << "NumberOfComponents=\"" << num_dim << "\" ";
  extrasStr << "format=\"" << "ascii" << "\"";
  Atz_XML_Package::writeTagStart(fid, nameTag.str().c_str(),
                                 extrasStr.str().c_str());

  // write vector data here
  for (int I = 0; I < numTotalMeshPts; I++) {
    for (int d = 0; d < num_dim; d++) {
      fid << vec_array[d][I] << " ";
    }
  }
  fid << endl;

  Atz_XML_Package::writeTagEnd(fid,"DataArray");

  nameTag.str("PointData");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  nameTag.str("CellData");
  Atz_XML_Package::writeTagStart(fid,nameTag.str().c_str());
  nameTag.str("CellData");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  nameTag.str("Piece");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  nameTag.str("RectilinearGrid");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  nameTag.str("VTKFile");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  fid.close();

  return 0;

}

int Atz_VTK_Helper::writeScalarFieldVTK_XML_File(const char   *filename,
                                                 int           num_dim,
                                                 int          *numMeshPtsPerDir,
                                                 double       *meshCenterX0,
                                                 double        meshDeltaX,
                                                 const char   *scalar_name,
                                                 double       *scalar_array) {

  const char *error_str_func = "writeScalarFieldVTK_XML_File()";

  stringstream extrasStr;
  stringstream nameTag;

  ofstream fid;

  // do some precalculations
  int N  = numMeshPtsPerDir[0]*numMeshPtsPerDir[1]*numMeshPtsPerDir[2];
  int Nw = numMeshPtsPerDir[0]*numMeshPtsPerDir[1]*(numMeshPtsPerDir[2]-1);

  int numTotalMeshPts = N;

  int J[3];
  int index;

  double deltaX;
  double X[3];

  // open the file
  fid.open(filename);

  if (!fid.is_open()) {  // if file is not open
    stringstream message;
    message << "Could not open file to write error ocurred." << endl;
    message << "  filename = " << filename << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  Atz_XML_Package::writeXMLHeader(fid);

  nameTag.str("VTKFile");
  extrasStr.str("");
  extrasStr << "type=\""       << "RectilinearGrid" << "\" ";
  extrasStr << "version=\""    << "0.1" << "\" ";
  extrasStr << "byte_order=\"" << "LittleEndian" << "\"";
  Atz_XML_Package::writeTagStart(fid,
                                 nameTag.str().c_str(),
                                 extrasStr.str().c_str());

  nameTag.str("RectilinearGrid");
  extrasStr.str("");
  extrasStr << "WholeExtent=\"" << "0 " << (numMeshPtsPerDir[0] - 1) << " "
                                << "0 " << (numMeshPtsPerDir[1] - 1) << " "
                                << "0 " << (numMeshPtsPerDir[2] - 1) << "\" ";
  extrasStr << "GhostLevel=\""  << "0" << "\" ";
  Atz_XML_Package::writeTagStart(fid,
                                 nameTag.str().c_str(),
                                 extrasStr.str().c_str());

  nameTag.str("Piece");
  extrasStr.str("");
  extrasStr << "Extent=\"" << "0 " << (numMeshPtsPerDir[0] - 1) << " "
                           << "0 " << (numMeshPtsPerDir[1] - 1) << " "
                           << "0 " << (numMeshPtsPerDir[2] - 1) << "\" ";
  Atz_XML_Package::writeTagStart(fid,
                                 nameTag.str().c_str(),
                                 extrasStr.str().c_str());

  nameTag.str("Coordinates");
  extrasStr.str("");
  extrasStr << "Scalars=\"" << scalar_name << "\"";
  Atz_XML_Package::writeTagStart(fid,
                                 nameTag.str().c_str(),
                                 extrasStr.str().c_str());

  // points along the x-axis, y-axis, z-axis
  for (int d = 0; d < num_dim; d++) {
    double L       = numMeshPtsPerDir[d]*meshDeltaX;
    double baseX_d = (meshCenterX0[d] - 0.5*L) + 0.5*meshDeltaX;
    //double baseX_d = (meshCenterX0[d] - 0.5*L);
    double X_d;
    double RangeMin = baseX_d;
    double RangeMax = baseX_d + (numMeshPtsPerDir[d] - 1)*meshDeltaX;


    nameTag.str("DataArray");
    extrasStr.str("");
    extrasStr << "type=\"" << "Float32" << "\" ";
    switch (d) {
    case 0:
      extrasStr << "Name=\"" << "x" << "\" ";
      break;
    case 1:
      extrasStr << "Name=\"" << "y" << "\" ";
      break;
    case 2:
      extrasStr << "Name=\"" << "z" << "\" ";
      break;
    default:
      extrasStr << "Name=\"" << "unknown-axis" << "\" ";
      break;
    } // end switch
    extrasStr << "format=\"" << "ascii" << "\" ";
    extrasStr << "RangeMin=\"" << RangeMin << "\" ";
    extrasStr << "RangeMax=\"" << RangeMax << "\"";
    Atz_XML_Package::writeTagStart(fid, nameTag.str().c_str(),
                                   extrasStr.str().c_str());

    // write the d-axis coordinates
    for (int m = 0; m < numMeshPtsPerDir[d]; m++) {
      X_d = baseX_d + m*meshDeltaX;
      fid << X_d << " ";
    }
    fid << endl;

    nameTag.str("DataArray");
    Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  } // end of the d-loop over different directions

  nameTag.str("Coordinates");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  nameTag.str("PointData");
  Atz_XML_Package::writeTagStart(fid,nameTag.str().c_str());
  // consider now each output array and if the flag is set

  nameTag.str("DataArray");
  extrasStr.str("");
  extrasStr << "type=\"" << "Float32" << "\" ";
  extrasStr << "Name=\"" << scalar_name << "\" ";
  extrasStr << "NumberOfComponents=\"" << 1 << "\" ";
  extrasStr << "format=\"" << "ascii" << "\"";
  Atz_XML_Package::writeTagStart(fid, nameTag.str().c_str(),
                                 extrasStr.str().c_str());

  // write vector data here
  for (int I = 0; I < numTotalMeshPts; I++) {
    fid << scalar_array[I] << " ";
  }
  fid << endl;

  Atz_XML_Package::writeTagEnd(fid,"DataArray");

  nameTag.str("PointData");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  nameTag.str("CellData");
  Atz_XML_Package::writeTagStart(fid,nameTag.str().c_str());
  nameTag.str("CellData");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  nameTag.str("Piece");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  nameTag.str("RectilinearGrid");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  nameTag.str("VTKFile");
  Atz_XML_Package::writeTagEnd(fid,nameTag.str().c_str());

  fid.close();

  return 0;

}








