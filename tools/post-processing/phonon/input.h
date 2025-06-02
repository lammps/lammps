#ifndef INPUT_H
#define INPUT_H

#include <cstdio>

class UserInput {
public:
   UserInput(int);
   ~UserInput();

   void read_stdin(char *);

private:
   FILE *fp;

};
#endif
