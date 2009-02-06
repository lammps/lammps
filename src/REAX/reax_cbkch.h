/////:EOH~
extern struct {
double ch[ReaxParams::nat];
} FORTRAN(cbkch,CBKCH);

extern struct {
  double chi[ReaxParams::nsort];
  double eta[ReaxParams::nsort];
  double gam[ReaxParams::nsort];
} FORTRAN(cbkchb,CBKCHB);
