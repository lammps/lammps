/* Numerical Recipes heap sort, modified to be like C library qsort */

/* accepts arbitrary objects to be sorted, user provides compare function
   this routine sorts a F77-style array indexed from 1-n
   thus you MUST call it from C with hpsort(ra-1, ) to offset C ptr by -1
   also added multiply-array-index-by-size to enable arbitrary objects
*/

void hpsort(char *ra, int n, int size,
	    int (*cmp)(const void *, const void *))
{
  unsigned long i,ir,j,l;
  char *rra;

  if (n < 2) return;
  l = (n >> 1)+1;
  ir = n;
  rra = (char *) malloc(size);
  for (;;) {
    if (l > 1) {
      memcpy(rra,&ra[(--l)*size],size);
    } else {
      memcpy(rra,&ra[ir*size],size);
      memcpy(&ra[ir*size],&ra[1*size],size);
      if (--ir == 1) {
	memcpy(&ra[1*size],rra,size);
	break;
      }
    }
    i = l;
    j = l+l;
    while (j <= ir) {
      if (j < ir && (*cmp)(&ra[j*size],&ra[(j+1)*size]) < 0) j++;
      if ((*cmp)(rra,&ra[j*size]) < 0) {
	memcpy(&ra[i*size],&ra[j*size],size);
	i = j;
	j <<= 1;
      } else j = ir+1;
    }
    memcpy(&ra[i*size],rra,size);
  }
  free(rra);
}
