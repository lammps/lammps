#ifndef PTM_NORMALIZE_VERTICES_H
#define PTM_NORMALIZE_VERTICES_H

namespace ptm {

void subtract_barycentre(int num, double (*points)[3], double (*normalized)[3]);
double normalize_vertices(int num, double (*points)[3], double (*normalized)[3]);

}

#endif

