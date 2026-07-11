#ifndef GREEN_H
#define GREEN_H

class Green{
public:
   Green(const int, const int, const int, const double, const double,
         const int, const double, double **, const int, double **);
   ~Green();

private:
   void Lanczos();
   void Recursion();
   void recursion();
 
   double **ldos;
 
   int natom, iatom, sysdim, nit, nw, ndim;
   double dw, wmin, wmax, epson;
   double **alpha, **beta, **H; 
   class Memory *memory;
};
#endif
