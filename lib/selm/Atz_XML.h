/*
 * Atz_XML.h
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#ifndef ATZ_XML_H_
#define ATZ_XML_H_

#include "mpi.h"
#include <string>
#include <map>

namespace Atz_XML {

  using namespace std;

  typedef map<string, string>   AttributesType;
  typedef pair<string, string>  AttributePairType;

  typedef map<string, void *>   ExtraDataType;
  typedef pair<string, void *>   ExtraDataPairType;

  static const int TAG_TYPE_NULL      = 0;
  static const int TAG_TYPE_START     = 1;
  static const int TAG_TYPE_END       = 2;
  static const int TAG_TYPE_EMPTY     = 3;
  static const int TAG_TYPE_QUESTMARK = 4;
  static const int TAG_TYPE_COMMENT   = 5;

}

#endif /* ATZ_XML_H_ */
