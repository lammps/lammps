/* ------------------------------------------------------------------------
 Provide some explanation here
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   This program is a courtesy of Aleksei Ivanov (Univ. of Iceland)
   Contributing authors: Aleksei Ivanov (Univ. of Iceland), 
   			 Julien Tranchida (SNL)
------------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

// calculate cubic coefficients

void count_coefficient(double *V, double *F, double *R, double *a, double *b,
    double *c,double *d,int M){
  /* R = square of distance between images*/
  /* V = energy of images */
  /* F = projection of real real forces along the path? */
  int i;
  for(i = 0; i < M ; i++) {
    a[i] = ( -2.0*(V[i+1]-V[i])/R[i] - F[i] - F[i+1] ) / (R[i]*R[i]);
    //a[i] = ( -2.0*(V[i+1]-V[i])/sqrt(R[i]) - F[i] - F[i+1] ) / R[i];
    b[i] = ( 3.0*(V[i+1]-V[i]) + (2.0*F[i]+F[i+1])*R[i] ) / (R[i]*R[i]);
    //b[i] = ( 3.0*(V[i+1]-V[i]) + (2.0*F[i]+F[i+1])*sqrt(R[i]) ) / R[i];
    c[i] = -F[i];
    d[i] = V[i];
  }
}

// cubic spline

double spline(double a,double b,double c,double d,double x) {
  return (a*x*x*x + b*x*x + c*x + d);
}

int main() {
  int M=0; 			// M+1 = number of images
  double *fmdottan;		// projection of real forces on tangent path
  double *coords;		// initial value of reaction coords
  double *V; 			// energy of images
  double *dist;			// square of the distance between images
  double *a, *b, *c, *d ;	// coefficients of cubic functions 
  double x;			// reaction coordinate
  double A,B; 			// additional variables for rnd
  double length = 0.0;
  int i,j;
  FILE *data;
  
  printf("Enter M = number of images - 1 \n");
  scanf("%d",&M);

  // allocating and initializing tables
  
  a = calloc(M,sizeof(double));			// cubic coefficients
  b = calloc(M,sizeof(double));
  c = calloc(M,sizeof(double));
  d = calloc(M,sizeof(double));
  V = calloc((M+1),sizeof(double)); 		// energies
  coords = calloc((M+1),sizeof(double));	// reaction coordinates 
  fmdottan = calloc((M+1),sizeof(double));	// fm dot tangent
  dist = calloc(M+1,sizeof(double));		// distance between images
  
  // reading input file

  if((data=fopen("neb_init.dat","r")) == NULL) {
    printf("Incorrect input file name.");
    return 0;
  }

  for(j=0; j < M+1; j++) {
    fscanf(data,"%lf\t%lf\t%lf\t%lf\n",&coords[j],&V[j],&fmdottan[j],&dist[j]);
    length += dist[j];
    printf("%lf %lf %lf %lf\n",coords[j],V[j],fmdottan[j],dist[j]);
  }

  if( (fclose(data)) == 0) {
  	printf("Data stored, input file closed.\n ");
  }

  // calculate value of coefficients
  
  count_coefficient(V,fmdottan,dist,a,b,c,d,M);
 
  // plot result of the interpolation

  if( ( data=fopen("interpolation_result.dat","w") )== NULL) {
    printf("Interpolation file could not be open.");
    return 0;
  }
 
  A = B = 0.0;
  for(i = 0; i < M ; i++) {
    B += dist[i];
    printf("%13le\n",B);
    for(j = 0; j <= 1000; j++) {
      x = dist[i]*1.0e-3*j;
      fprintf(data,"%13lf\t%13le\n",(x+A)/length,spline(a[i],b[i],c[i],d[i],x));
    }
    A += dist[i];
  }

  return 0;
}
