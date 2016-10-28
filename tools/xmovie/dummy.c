/* program to generate dummy data set */
/* Mike Uttormark 7/16/92 */
/* Sandia Nat'l Labs 1421 */
/* On leave from University of Wisconsin--Madison */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define RANGE   5.0
#define STEPS	2000
#define PI	3.14159
#define S	50

#ifdef MISSINGDEFS
int printf(char *, ...);
#endif

int main(int argc, char **argv);

int main(int argc, char **argv)
{
	int	i, j, step, npart;
	int	nsteps;

	if (argc < 2)
		nsteps = STEPS;
	else
		nsteps = atoi(argv[1]);

	if (argc < 3)
		npart = 1;
	else
		npart = atoi(argv[2]);

	for(step = 0; step < nsteps; step++) {

		printf("ITEM: TIME\n%g\n", (float) step);

		if (step == 0) {
			printf("ITEM: BOUNDS\n");
			for (i = 3; i; i--)
				printf("%g, %g\n", -RANGE, RANGE);
		}
		printf("ITEM: POSITIONS\n");

		for( j = 0; j < npart; j++)
			printf("%i, %i, %g, %g, %g\n",
				j+1,
				1 + j%4,
				RANGE*sin(2*PI*(step+j)/S),
				RANGE*cos(0.25*PI*(step+j)/S),
				RANGE*sin(5*PI*(step+j)/S));

		if (npart > 1) {
			printf("ITEM: BONDS\n");
			for ( j = 0; j < npart-1; j++)
				printf("%i %i %i\n", 1 + j%4, j+1, j+2);
		}
	}
}


