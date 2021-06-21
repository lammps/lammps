/* ----------------------------------------------------------------------

 Parser
 
 Paul J. Atzberger
 http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_PARSER1_H
#define SELM_PARSER1_H

namespace LAMMPS_NS {

class SELM_Parser1 {

 public:

  /* ================ Constants ================= */
  static const int   MAX_NUM_PARAMS        = 100;
  static const int   MAX_NUM_FIDS          = 10;

  static const int   PARAMTYPE_NULL        = 0;
  static const int   PARAMTYPE_INT         = 1;
  static const int   PARAMTYPE_DOUBLE      = 2;
  static const int   PARAMTYPE_STRING      = 3;
  static const int   PARAMTYPE_CHAR        = 4;
  static const int   PARAMTYPE_INT_LIST    = 5;
  static const int   PARAMTYPE_DOUBLE_LIST = 6;
  static const int   PARAMTYPE_STRING_LIST = 7;  /* Note list must be terminated by \0*/

  static const char* KEYWORD_INCLUDE;

  /* ================ Data structure type definitions ================= */
  typedef struct paramDescrType {
    char paramName[100];
    int  paramType;
    int  paramSetFlag;  /* whethor parameter found in file or not */
    void *paramVar;
    void *paramExtras;
  } paramDescrType;

  typedef struct paramSpecificationType {
    int numParams;
    paramDescrType *paramDescrList;
  } paramSpecificationType;


  /* ================ Function prototypes ================= */
  SELM_Parser1();
  virtual ~SELM_Parser1();

  void parseParameters(const char *filename, paramSpecificationType *paramSpecification);
  void printParameters(paramSpecificationType *paramSpecification);
  int  areAllParametersSet(paramSpecificationType *paramSpecification);
  void printUnsetParameters(paramSpecificationType *paramSpecification);

  void packageError(int errorCode, void *extras);

};

}

#endif
