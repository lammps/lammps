// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file common.hh
 * \brief Header file for the small helper functions. */

#ifndef VOROPP_COMMON_HH
#define VOROPP_COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "config.hh"

namespace voro {

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
inline void voro_fatal_error(const char *p,int status) {
	fprintf(stderr,"voro++: %s\n",p);
	exit(status);
}

/** \brief Prints a vector of positions.
 *
 * Prints a vector of positions as bracketed triplets.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
inline void voro_print_positions(std::vector<double> &v,FILE *fp=stdout) {
	if(v.size()>0) {
		fprintf(fp,"(%g,%g,%g)",v[0],v[1],v[2]);
		for(int k=3;(unsigned int) k<v.size();k+=3) {
			fprintf(fp," (%g,%g,%g)",v[k],v[k+1],v[k+2]);
		}
	}
}

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation
 * was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
inline FILE* safe_fopen(const char *filename,const char *mode) {
	FILE *fp=fopen(filename,mode);
	if(fp==NULL) {
		fprintf(stderr,"voro++: Unable to open file '%s'\n",filename);
		exit(VOROPP_FILE_ERROR);
	}
	return fp;
}

void voro_print_vector(std::vector<int> &v,FILE *fp=stdout);
void voro_print_vector(std::vector<double> &v,FILE *fp=stdout);
void voro_print_face_vertices(std::vector<int> &v,FILE *fp=stdout);

}

#endif
