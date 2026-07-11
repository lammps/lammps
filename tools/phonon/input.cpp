#include "input.h"

#include "global.h"

/* -------------------------------------------------------------------
 * Constructor. If flag = 1, output user inputs as script.inp
 * ---------------------------------------------------------------- */ 
UserInput::UserInput(int flag)
{
   fp = nullptr;
   if (flag) fp = fopen("script.inp", "w");
}

/* -------------------------------------------------------------------
 * Deconstructor. Output user inputs as required and clear workspace.
 * ---------------------------------------------------------------- */ 
UserInput::~UserInput()
{
   if (fp) fclose(fp);
   fp = nullptr;
}

/* -------------------------------------------------------------------
 * Read stdin and keep a record of it.
 * ---------------------------------------------------------------- */ 
void UserInput::read_stdin(char *str)
{
   fgets(str, MAXLINE, stdin);
   if (fp) fprintf(fp, "%s", str);
}
/* ---------------------------------------------------------------- */
