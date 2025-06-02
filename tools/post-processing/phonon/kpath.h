#ifndef UseSPG
#define KPATH_H
#endif
#ifndef KPATH_H
#define KPATH_H

class kPath{
public:
  kPath(class DynMat *, class QNodes *);
  ~kPath();

  void kpath();
  void show_path();
  void show_info();

private:
  class Memory *memory;
  class DynMat *dynmat;
  class QNodes *q;
  char symbol[11];
  int spgnum, sysdim, num_atom, *attyp;
  double latvec[3][3], **atpos;

};
#endif
