/////:EOH~
extern struct {
  int nvl1[ReaxParams::nneighmax * ReaxParams::nat];
  int nvl2[ReaxParams::nneighmax * ReaxParams::nat];
  int nvpair;
  int nvlself;
} FORTRAN(cbkpairs,CBKPAIRS);

extern struct {
  int nvlbo[ReaxParams::nneighmax * ReaxParams::nat];
} FORTRAN(cbknvlbo,CBKNVLBO);

extern struct {
  int nvlown[ReaxParams::nneighmax * ReaxParams::nat];
} FORTRAN(cbknvlown,CBKNVLOWN);
