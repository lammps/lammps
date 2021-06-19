/* ----------------------------------------------------------------------

   Stochastic Eulerian Lagrangian Methods

   Paul J. Atzberger
   http://atzberger.org/
 

------------------------------------------------------------------------- */

#include "Atz_XML_Package.h"

using namespace std; /* ensures standard template library in namespace */

/* ---------------------------------------------------------------------- */
/* Constant valued strings                                                */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

/* initialize static variables */

Atz_XML_Package::Atz_XML_Package() {

  const char *error_str_code = "Atz_XML_Package.cpp";
  const char *error_str_func = "Atz_XML_Package()";

}

Atz_XML_Package::~Atz_XML_Package() {

}


void Atz_XML_Package::writeXMLHeader(FILE *fid) {
  stringstream output;

  output << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;

  /* write the output to disk */
  fprintf(fid,"%s",output.str().c_str());

}

void Atz_XML_Package::writeXMLHeader(ofstream &fid) {
  fid << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
}

void Atz_XML_Package::writeTagStart(FILE *fid, const char *tagName) {
  writeTagStart(fid,tagName,"");
}

void Atz_XML_Package::writeTagStart(FILE *fid, string tagName) {
  writeTagStart(fid, tagName.c_str());
}

void Atz_XML_Package::writeTagStart(ofstream &fid, const char *tagName) {
  writeTagStart(fid,tagName,"");
}

void Atz_XML_Package::writeTagStart(ofstream &fid, string tagName) {
  writeTagStart(fid, tagName.c_str());
}

void Atz_XML_Package::writeTagStart(ofstream &fid, string tagName, string extras) {
  writeTagStart(fid, tagName.c_str(), extras.c_str());
}

void Atz_XML_Package::writeTagStart(FILE *fid, const char *tagName, const char *extras) {

  stringstream output;

  if (extras[0] == 0) {
    output << "<" << tagName << ">" << endl;
  } else {
    output << "<" << tagName << " " << extras << ">" << endl;
  }

  /* write the output to disk */
  fprintf(fid,"%s",output.str().c_str());

}

void Atz_XML_Package::writeTagStart(ofstream &fid, const char *tagName, const char *extras) {

  if (extras[0] == 0) {
    fid << "<" << tagName << ">" << endl;
  } else {
    fid << "<" << tagName << " " << extras << ">" << endl;
  }

}

void Atz_XML_Package::writeTagEnd(FILE *fid, string tagName) {
  writeTagEnd(fid, tagName.c_str());
}

void Atz_XML_Package::writeTagEnd(ofstream &fid, string tagName) {
  writeTagEnd(fid, tagName.c_str());
}

void Atz_XML_Package::writeTagEnd(FILE *fid, const char *tagName) {

  stringstream output;

  output << "</" << tagName << ">" << endl;

  /* write the output to disk */
  fprintf(fid,"%s",output.str().c_str());

}

void Atz_XML_Package::writeTagEnd(ofstream &fid, const char *tagName) {
  fid << "</" << tagName << ">" << endl;
}

void Atz_XML_Package::writeTagValueDouble(FILE *fid, const char *tagName, double value) {
  writeTagValueDoubleArray(fid, tagName, 1, &value);
}

void Atz_XML_Package::writeTagValueDouble(ofstream &fid, const char *tagName, double value) {
  writeTagValueDoubleArray(fid, tagName, 1, &value);
}

void Atz_XML_Package::writeTagValueDoubleArray(FILE *fid, const char *tagName, int N, double *values) {

  stringstream output;

  output << "<" << tagName << " value="
         << "\"";
  for(int i = 0; i < N; i++) {
    output << values[i];
    if (i != N - 1)
      output << " ";  // insert space (expect last entry)
  }
  output << "\"" << "/>" << endl;

  /* write the output to disk */
  fprintf(fid,"%s",output.str().c_str());

}

void Atz_XML_Package::writeTagValueDoubleArray(ofstream &fid, const char *tagName, int N, double *values) {

  fid << "<" << tagName << " value="
         << "\"";
  for(int i = 0; i < N; i++) {
    fid << values[i];
    if (i != N - 1)
      fid << " ";  // insert space (expect last entry)
  }
  fid << "\"" << "/>" << endl;

}



void Atz_XML_Package::packageError(const char *error_str_code, const char *error_str_func, stringstream &message) {
  packageError(error_str_code, error_str_func, message.str().c_str());
}


void Atz_XML_Package::packageError(const char *error_str_code, const char *error_str_func, string &message) {
  packageError(error_str_code, error_str_func, message.c_str());
}


void Atz_XML_Package::packageError(const char *error_str_code, const char *error_str_func, const char *message) {
  SELM_Package::packageError(error_str_code, error_str_func, message);
}
