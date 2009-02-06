/////:EOH~
extern struct {
  double d[3*ReaxParams::nat]; double estrain[ReaxParams::nat];
} FORTRAN(cbkd,CBKD);

extern struct {
  double atomvirial[6*ReaxParams::nat];
  double virial[6]; 
  int Lvirial;
  int Latomvirial;
} FORTRAN(cbkvirial,CBKVIRIAL);

