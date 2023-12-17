/****************************
*
* This function first allocates memory to the forcefield item
* structures and then reads parameters from the forcefield file into the
* allocated memory
*
*/

#include "msi2lmp.h"
#include "Forcefield.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#if defined(_WIN32)
#define strdup(x) _strdup(x)
#endif

static int blank_line(char *line)
{
  while (*line != '\0') {
    if (isalnum((int) *line)) return 0;
    ++line;
  }
  return 1;
}

static unsigned char string_match(const char *,const char *);

void ClearFrcItem(struct FrcFieldItem *item)
{
    free(item->data);
}

const char *SearchAndCheck(const char *keyword)
{
  char *status;
  int got_it = 0;
  char line[MAX_LINE_LENGTH] = "empty";

  rewind(FrcF);
  while (got_it == 0) {
    status = fgets( line, MAX_LINE_LENGTH, FrcF );
    if (status == NULL) {
      fprintf(stderr," Unable to find keyword '%s'\n",keyword);
      fprintf(stderr," Check consistency of forcefield name and class \n");
      fprintf(stderr," Exiting....\n");
      exit(1);
    }
    if (has_utf8(line)) utf8_subst(line);
    if (line[0] == '@') {
      if (string_match(strtok(line+1," '\t\n\r\f("),keyword)) {
        got_it = 1;
        status = strtok(NULL," '\t\n\r\f(");
        if (status != NULL)
          return strdup(status);
      }
    }
  }
  return strdup("(unknown)");
}

void SearchAndFill(struct FrcFieldItem *item)
{
  int i,j;  /* counters */
  int got_it = 0;
  int ctr = 0;
  long file_pos;
  char line[MAX_LINE_LENGTH] = "empty";
  char *charptr,*status;

  /***********************ALLOCATE MEMORY FOR STRUCTURE ********************/

  /* Read and discard lines until keyword is found */

  rewind(FrcF);
  while (got_it == 0) {
    status = fgets( line, MAX_LINE_LENGTH, FrcF );
    if (status == NULL) {
      fprintf(stderr," Unable to find keyword '%s'\n",item->keyword);
      fprintf(stderr," Check consistency of forcefield name and class \n");
      fprintf(stderr," Exiting....\n");
      exit(1);
    }
    if (has_utf8(line)) utf8_subst(line);
    if (line[0] == '#') {
      if (string_match(strtok(line," '\t\r\n("),item->keyword)) got_it = 1;
    }
    /*     if (strncmp(line, item->keyword,strlen(item->keyword))==0) got_it = 1; */
  }

  file_pos = ftell(FrcF);
  if (file_pos < 0) {
    fprintf(stderr, "Could not obtain file stream position: %s\n", strerror(errno));
    exit(2);
  }

  /* Count the number of lines until next item is found */

  while (strncmp(fgets(line,MAX_LINE_LENGTH,FrcF), "#", 1) != 0 )
    ctr++;

  /* Allocate the memory using calloc */

  item->data = (struct FrcFieldData *)calloc(ctr, sizeof(struct FrcFieldData));

  if (item->data == NULL) {
    fprintf(stderr,"Could not allocate memory to %s\n", item->keyword);
    exit(2);
  }

  /********************FILL PARAMETERS AND EQUIVALENCES ********************/

  /* Read lines until keyword is found */

  if (fseek(FrcF,file_pos,SEEK_SET) < 0) {
    fprintf(stderr, "Resetting file stream failed: ", strerror(errno));
    exit(2);
  }
  strcpy(line,"empty");

  /* Read lines until data starts (when !--- is found) */

  ctr = 0;
  while ( strncmp(line,"!---", 4) != 0 ) {
    fgets(line, MAX_LINE_LENGTH, FrcF);
    if (has_utf8(line)) utf8_subst(line);
  }

  /* Get first line of data that isn't commented out */

  fgets(line, MAX_LINE_LENGTH, FrcF);
  if (has_utf8(line)) utf8_subst(line);
  while (strncmp(line,"!",1) == 0) {
    fgets( line, MAX_LINE_LENGTH, FrcF);
    if (has_utf8(line)) utf8_subst(line);
  }

  /* Read data into structure */

  while( strncmp( line, "#", 1 ) != 0 ) {

    float version;
    int reference,replace;
    char atom_types[5][5];
    double parameters[8];

    /* version number and reference number */

    version = atof(strtok(line, WHITESPACE));
    reference = atoi(strtok(NULL, WHITESPACE));

    /* equivalences */

    for(i = 0; i < item->number_of_members; i++ ) {
      charptr = strtok(NULL, WHITESPACE);
      if (strlen(charptr) > 4) {
        fprintf(stderr,"Warning: type name overflow for '%s'. "
                "Truncating to 4 characters.\n",charptr);
      }
      sscanf(charptr,"%4s",atom_types[i]);
    }

    /* parameters -- Because of symmetrical terms, bonang, angtor, and
       endbontor have to be treated carefully */

    for( i = 0; i < item->number_of_parameters; i++ ) {
      charptr = strtok(NULL, WHITESPACE);
      if(charptr == NULL) {
        for ( j = i; j < item->number_of_parameters; j++ )
          parameters[j] = parameters[j-i];
        break;
      } else {
        parameters[i] = atof(charptr);
      }
    }
    /* Search for matching sets of atom types.
       If found and the version number is greater, substitute
       the current set of parameters in place of the found set.
       Otherwise, add the current set of parameters to the
       list.
    */
    replace = ctr;
    for (j=0; j < ctr; j++) {

      int k=0;
      int match = 1;
      while (match && (k < item->number_of_members)) {
        if (strncmp(item->data[j].ff_types[k],atom_types[k],5) == 0)
          k++;
        else
          match = 0;
      }
      if (match == 1) {
        replace = j;
        break;
      }
    }
    if (replace != ctr) {
      if (version > item->data[replace].ver) {

        if (pflag > 1) {
          fprintf(stderr," Using higher version of parameters for");
          fprintf(stderr," %s  ",item->keyword);
          for (i=0; i < item->number_of_members; i++)
            fprintf(stderr,"%s ",atom_types[i]);
          fprintf(stderr," version %3.2f\n",version);
        }

        item->data[replace].ver = version;
        item->data[replace].ref = reference;
        for (i=0; i < item->number_of_members; i++) {
          strncpy(item->data[replace].ff_types[i],atom_types[i],5);
        }
        for (i=0; i < item->number_of_parameters; i++) {
          item->data[replace].ff_param[i] = parameters[i];
        }
      } else {
        if (pflag > 1) {
          fprintf(stderr," Using higher version of parameters for");
          fprintf(stderr," %s  ",item->keyword);
          for (i=0; i < item->number_of_members; i++)
            fprintf(stderr,"%s ",item->data[replace].ff_types[i]);
          fprintf(stderr," version %3.2f\n",item->data[replace].ver);
        }
      }
    } else {
      item->data[ctr].ver = version;
      item->data[ctr].ref = reference;
      for (i=0; i < item->number_of_members; i++) {
        strncpy(item->data[ctr].ff_types[i],atom_types[i],5);
      }
      for (i=0; i < item->number_of_parameters; i++) {
        item->data[ctr].ff_param[i] = parameters[i];
      }
      ctr++;
    }
    fgets( line, MAX_LINE_LENGTH, FrcF);
    if (has_utf8(line)) utf8_subst(line);

    /*if blank line encountered, get next */
    while((blank_line(line)) || (strncmp(line,"!",1) == 0)) {
      status = fgets( line, MAX_LINE_LENGTH, FrcF);
      if (status == NULL) break;
      if (has_utf8(line)) utf8_subst(line);
    }
  }
  item->entries = ctr;

  /*Debugging
    fprintf(stderr,"\n%s\n", item->keyword);
    for(i=0;i<ctr;i++) {
    for(j=0;j<item->number_of_members;j++)
    fprintf(stderr,"%3s ", item->data[i].ff_equiv[j]);
    fprintf(stderr,"     ");
    for(j=0;j<item->number_of_parameters;j++)
    fprintf(stderr,"%10.5f ",item->data[i].ff_param[j]);
    fprintf(stderr,"\n");
    }
  */
}

unsigned char string_match(const char *string1,const char *string2)
{
  int len1,len2;

  len1 = strlen(string1);
  len2 = strlen(string2);

  if (len1 != len2) {
    return 0;
  } else {
    if (strncmp(string1,string2,len1) == 0) {
      return 1;
    } else {
      return 0;
    }
  }
}
