/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 
------------------------------------------------------------------------- */

#include <cstdlib>
#include "stdio.h"
#include "string.h"
#include "SELM_Parser1.h"

#include "error.h"

#include <limits.h>

#include <malloc.h>

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */
/* Constant valued strings                                                */
/* ---------------------------------------------------------------------- */

/* define string constants */
const char* SELM_Parser1::KEYWORD_INCLUDE = "@include";

/* integer constants were defined within the class */

/* ---------------------------------------------------------------------- */

SELM_Parser1::SELM_Parser1() {

  const char *error_str_code = "SELM_Parser1.cpp";
  const char *error_str_func = "SELM_Parser1()";

}

SELM_Parser1::~SELM_Parser1() {

}

void SELM_Parser1::parseParameters(const char *filename,
                                   paramSpecificationType *paramSpecification) {

  const char *error_str_code = "SELM_Parser1.cpp";
  const char *error_str_func = "parseParameters()";

  FILE *fid;
  FILE *fid_list[MAX_NUM_FIDS];
  int   fidListIndex = 0;
  int   quitFlag = 0;
  char  inputString[10000] = "";
  int   stringNumber = 0;
  int   scanFlag;

  int   verboseFlag = 0;

  int   listFlag  = 0;
  int   listIndex = 0;
  int   listType  = 0;

  int   paramMatchFlag = 0;

  int   i,j,k;

  int   c;

  int   matchFlag;

  int   paramIndex;
  int   listParamIndex;

  int    tmpInt;
  double tmpDouble;

  /* Open the file and record it in the fid list */
  fid = fopen(filename,"r");
  fid_list[fidListIndex] = fid;
  fidListIndex++;

  if (fid == 0) {
    printf("Invalid Filename %s \n",filename);
    exit(1);
  }

  while (!quitFlag) {
    scanFlag = fscanf(fid,"%s",&inputString);
    if (scanFlag != EOF) {
      if (verboseFlag)
        printf("[%d] %d: %s \n",scanFlag,stringNumber,inputString);
      stringNumber++;

      if (inputString[0] == '#') { /* check if line is a comment */
        if (verboseFlag)
          printf("Line is a comment \n");
        c = fgetc(fid);
        while ((c != EOF) && (c != '\n')) {/* loops to read characters
                          after # */
          c = fgetc(fid);
          if (verboseFlag)
            printf(" Comment Char: %1s \n",&c);
        }

        if (c == EOF) {
          quitFlag = 1;
        }

      } else { /* Check string against parameter names */
        matchFlag = 0;
        if ((listType == PARAMTYPE_STRING_LIST) && (listFlag == 1)) {
          /* Skip the check if in string list mode */
        } else {
          /* Check against the keyword names */
          if (strcmp(KEYWORD_INCLUDE,
                     inputString) == 0) {
            if (verboseFlag) {
              printf("Detected Keyword %s \n",KEYWORD_INCLUDE);
            }

            scanFlag = fscanf(fid,
                              "%s",
                              inputString);

            if (verboseFlag)
              printf(" Attempting to include contents of file %s \n",
                     inputString);

            /* Open the file */
            fid = fopen(inputString,"r");
            if (fid == 0) {
              printf(" Error Opening file %s \n",inputString);
              exit(1);
            }

            /* Add fid to the file list */
            if (fidListIndex < MAX_NUM_FIDS) {
              fid_list[fidListIndex] = fid;
              fidListIndex++;
            } else {
              printf("Maximum number of @includes invoked, %d \n",
                     MAX_NUM_FIDS);
              printf("May have included its own file \n");
              exit(1);
            }

            /* Proceed using this file until EOF */
            /* Revert to last fid when reach EOF */
          }

          /* Check against the parameter names */
          for (k = 0; k < paramSpecification->numParams; k++) {
            if (strcmp(paramSpecification->paramDescrList[k].paramName,
                       inputString) == 0) {
              if (verboseFlag)
                printf("Detected Parameter %s \n",
                       paramSpecification->paramDescrList[k].paramName);

              switch(paramSpecification->paramDescrList[k].paramType) {

              case PARAMTYPE_INT:

                if (verboseFlag)
                  printf(" Parameter is of type int \n");

                scanFlag = fscanf(fid,
                                  "%s",
                                  inputString);

                if (verboseFlag)
                  printf(" Attempting to convert %s to int \n",inputString);

                /* old way of doing the parsing */
                /* *(int *)paramSpecification->paramDescrList[k].paramVar = atoi(inputString); */

                /* Use the ASCII conversion to double to allow for scientific notation and
                 * use type casting to obtain the specific number type */
                tmpDouble = strtod(inputString, NULL);
                if ((tmpDouble > INT_MAX) || (tmpDouble < INT_MIN)) { /* check value is within range */
                  printf("ERROR: %s, %s \n", error_str_code, error_str_func);
                  printf("Integer value is outside of range of the C int type. \n");
                  printf("Value specified = %s \n", inputString);
                  printf("INT_MAX = %d \n", INT_MAX);
                  printf("INT_MIN = %d \n", INT_MIN);
                  printf("paramName = %s \n",paramSpecification->paramDescrList[k].paramName);
                  packageError(1,0);
                }
                tmpInt    = (int)tmpDouble; /* convert the value to an integer */

                *(int *)paramSpecification->paramDescrList[k].paramVar = tmpInt; /* save the value */

                paramSpecification->paramDescrList[k].paramSetFlag = 1;
                matchFlag = 1;

              break;

              case PARAMTYPE_STRING:

                if (verboseFlag)
                  printf(" Parameter is of type string \n");

                scanFlag = fscanf(fid,
                                  "%s",
                                  paramSpecification
                                  ->paramDescrList[k].paramVar);

                paramSpecification->paramDescrList[k].paramSetFlag = 1;

                matchFlag = 1;

              break;

              case PARAMTYPE_DOUBLE:

                if (verboseFlag)
                  printf(" Parameter is of type double \n");

                scanFlag = fscanf(fid,
                                  "%s",
                                  inputString);

                if (verboseFlag)
                  printf(" Attempting to conver %s to double \n",inputString);


                /* Use the ASCII conversion to double to allow for scientific notation and
                 * use type casting to obtain the specific number type */
                tmpDouble = strtod(inputString, NULL);

                *(double *)paramSpecification->paramDescrList[k].paramVar = tmpDouble;  /* atof(inputString) */

                paramSpecification->paramDescrList[k].paramSetFlag = 1;

                matchFlag = 1;

              break;

              case PARAMTYPE_INT_LIST:

                if (verboseFlag) {
                  printf(" Parameter is of type int list \n");
                  printf(" Setting List flags \n");
                }

                listFlag       = 1;
                listParamIndex = k;
                listType       = PARAMTYPE_INT_LIST;
                listIndex      = 0;

                matchFlag      = 2;

                break;

              case PARAMTYPE_DOUBLE_LIST:

                if (verboseFlag) {
                  printf(" Parameter is of type double list \n");
                  printf(" Setting List flags \n");
                }

                listFlag       = 1;
                listParamIndex = k;
                listType       = PARAMTYPE_DOUBLE_LIST;
                listIndex      = 0;
                matchFlag      = 2;

              break;

              case PARAMTYPE_STRING_LIST:

                if (verboseFlag) {
                  printf(" Parameter is of type string list \n");
                  printf(" Setting List flags \n");
                }

                listFlag       = 1;
                listParamIndex = k;
                listType       = PARAMTYPE_STRING_LIST;
                listIndex      = 0;

                matchFlag      = 2;
              break;

              }
            }
          }

          /* Display an error message if the inputString can not be
             matched against a parameter name */
          /*
             if (matchFlag == 0) {
               printf("PARSE ERROR (paramParse.c): The string %s could not be ");
               printf("matched with a parameter name in parse list.\n",
               inputString);
               exit(1);
             }
           */

        }

        if (listFlag == 1) {

          /* For the string list only consider an end of the
             list if \0 string encountered.  This allows
             for keywords to appear in the string list */
          if (listType == PARAMTYPE_STRING_LIST) {

            if (strcmp(inputString,"\\0") != 0) {

              if (matchFlag == 1) {
                /* modify match to treat keyword as list element */
                matchFlag = 0;
                if (verboseFlag)
                  printf("Matched Keyword, but ignoring as string list \n");
              }

            } else {

              if (verboseFlag)
                printf("Encountered String list termination \\0 \n");

              paramSpecification->paramDescrList[listParamIndex].paramSetFlag
              = 1;

              matchFlag = 1; /* signal to end the list */

            }
          }


          /* For all lists consider if keyword matched or not in
             scan of the last inputString (string list modifies flag) */
          switch (matchFlag) {

          case 0: /* If no match then use word to build list */

            if (listType == PARAMTYPE_INT_LIST) {

              if (verboseFlag)
                printf(" Attempting to convert %s to int \n",inputString);

              /* Use the ASCII conversion to double to allow for scientific notation and
               * use type casting to obtain the specific number type */
              tmpDouble = strtod(inputString, NULL);
              if ((tmpDouble > INT_MAX) || (tmpDouble < INT_MIN)) { /* check value is within range */
                printf("ERROR: %s, %s \n", error_str_code, error_str_func);
                printf("Integer value is outside of range of the C int type. \n");
                printf("Value specified = %s \n", inputString);
                printf("INT_MAX = %d \n", INT_MAX);
                printf("INT_MIN = %d \n", INT_MIN);
                printf("paramName = %s \n",paramSpecification->paramDescrList[k].paramName);
                packageError(1,0);
              }
              tmpInt    = (int)tmpDouble; /* convert the value to an integer */

              *(((int *)paramSpecification->
                  paramDescrList[listParamIndex].paramVar)
                  + listIndex)
                  = tmpInt;

              listIndex++;

              *(int *)paramSpecification->
                  paramDescrList[listParamIndex].paramExtras
                  = listIndex; /* extra parameter holds number of elements */

              paramSpecification->paramDescrList[listParamIndex].paramSetFlag
              = 1;

            }

            if (listType == PARAMTYPE_DOUBLE_LIST) {

              if (verboseFlag)
                printf(" Attempting to convert %s to double \n",inputString);

              /* Use the ASCII conversion to double to allow for scientific notation and
               * use type casting to obtain the specific number type */
              tmpDouble = strtod(inputString, NULL);
              *(((double *)paramSpecification->
                  paramDescrList[listParamIndex].paramVar)
                  + listIndex)
                  = tmpDouble;

              listIndex++;

              *(int *)paramSpecification->
                  paramDescrList[listParamIndex].paramExtras
                  = listIndex; /* extra parameter holds number of elements */

              paramSpecification->paramDescrList[listParamIndex].paramSetFlag
              = 1;

            }

            if (listType == PARAMTYPE_STRING_LIST) {

              if (verboseFlag)
                printf(" Attempting to put %s in list  \n",inputString);

              strcpy(((char **)paramSpecification->
                  paramDescrList[listParamIndex].paramVar)[listIndex]
                                                           ,inputString);
              listIndex++;

              *(int *)paramSpecification->
                  paramDescrList[listParamIndex].paramExtras
                  = listIndex; /* extra parameter holds number of elements */

              paramSpecification->paramDescrList[listParamIndex].paramSetFlag
              = 1;
            }

            break;

          case 1: /* Found match to word, if building list then end it */
            listFlag = 0;
            break;

          case 2:
            /* Wait until next word read to begin building the list */
            break;

          }
        }
      }

    } else {
      fidListIndex--;

      if (fidListIndex == 0) {
        quitFlag = 1;
        if (verboseFlag)
          printf("End of Root File \n");
      } else {
        /* Close File */
        fclose(fid);

        /* Revert to already opened file */
        fid = fid_list[fidListIndex];
        if (verboseFlag) {
          printf("EOF for fidIndex %d reverting to fidIndex %d \n",
                 fidListIndex + 1,fidListIndex);
        }
      }

    }
  }

  /* Close Files */
  fclose(fid);

}


void SELM_Parser1::printParameters(paramSpecificationType *paramSpecification) {

  int j,k;

  for (k = 0; k < paramSpecification->numParams; k++) {
    switch(paramSpecification->paramDescrList[k].paramType) {
    case PARAMTYPE_INT:
      printf(" %s = %d \n",paramSpecification->paramDescrList[k].paramName,
             *(int *)paramSpecification->paramDescrList[k].paramVar);
      break;

    case PARAMTYPE_DOUBLE:
      printf(" %s = %g \n",paramSpecification->paramDescrList[k].paramName,
             *(double *)paramSpecification->paramDescrList[k].paramVar);
      break;

    case PARAMTYPE_STRING:
      printf(" %s = %s \n",paramSpecification->paramDescrList[k].paramName,
             (char *)paramSpecification->paramDescrList[k].paramVar);
      break;

    case PARAMTYPE_INT_LIST:
      printf(" %s.N = %d \n",paramSpecification->paramDescrList[k].paramName,
             *(int *)paramSpecification->paramDescrList[k].paramExtras);
      for (j = 0;
          j < *(int *)paramSpecification->paramDescrList[k].paramExtras;
          j++) {
        printf(" %s.list[%d] = %d \n",
               paramSpecification->paramDescrList[k].paramName,
               j,
               ((int *)(paramSpecification->paramDescrList[k].paramVar))[j]);
      }
      break;

    case PARAMTYPE_DOUBLE_LIST:
      printf(" %s.N = %d \n",paramSpecification->paramDescrList[k].paramName,
             *(int *)paramSpecification->paramDescrList[k].paramExtras);
      for (j = 0;
          j < *(int *)paramSpecification->paramDescrList[k].paramExtras;
          j++) {
        printf(" %s.list[%d] = %g \n",
               paramSpecification->paramDescrList[k].paramName,
               j,
               ((double *)(paramSpecification->paramDescrList[k].paramVar))[j]);
      }
      break;

    case PARAMTYPE_STRING_LIST:
      printf(" %s.N = %d \n",paramSpecification->paramDescrList[k].paramName,
             *(int *)paramSpecification->paramDescrList[k].paramExtras);
      for (j = 0;
          j < *(int *)paramSpecification->paramDescrList[k].paramExtras;
          j++) {
        printf(" %s.list[%d] = %s \n",
               paramSpecification->paramDescrList[k].paramName,
               j,
               ((char *)paramSpecification->paramDescrList[k].paramVar)[j]);
      }
      break;

    }
  }

}


int SELM_Parser1::areAllParametersSet(paramSpecificationType *paramSpecification) {

  int j,k;
  int areAllSet = 1;

  for (k = 0; k < paramSpecification->numParams; k++) {
    areAllSet = areAllSet
        && paramSpecification->paramDescrList[k].paramSetFlag;
  }

  return areAllSet;

}


void SELM_Parser1::printUnsetParameters(paramSpecificationType *paramSpecification) {

  int j,k;
  int areAllSet = 1;
  char typeString[100];

  for (k = 0; k < paramSpecification->numParams; k++) {
    if (paramSpecification->paramDescrList[k].paramSetFlag == 0) {
      switch(paramSpecification->paramDescrList[k].paramType) {
      case PARAMTYPE_INT:
        strcpy(typeString,"Integer");
        break;
      case PARAMTYPE_DOUBLE:
        strcpy(typeString,"Double");
        break;
      case PARAMTYPE_STRING:
        strcpy(typeString,"String");
        break;
      case PARAMTYPE_INT_LIST:
        strcpy(typeString,"Integer List");
        break;
      case PARAMTYPE_DOUBLE_LIST:
        strcpy(typeString,"Double List");
        break;
      case PARAMTYPE_STRING_LIST:
        strcpy(typeString,"String List (\0 termination tag)");
        break;
      }
      printf("%s : %s \n",
             paramSpecification->paramDescrList[k].paramName,
             typeString);
    }
  }

}


void SELM_Parser1::packageError(int code, void *extras) {

    exit(code);
}
