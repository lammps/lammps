#include "stdio.h"
#include "stdlib.h"
#include "dynmat.h"
#include "phonon.h"

using namespace std;

int main(int argc, char** argv)
{

  DynMat *dynmat = new DynMat(argc, argv);
  Phonon *phonon = new Phonon(dynmat);

  delete phonon;
  delete dynmat;

return 0;
}
