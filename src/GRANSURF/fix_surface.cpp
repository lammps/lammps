/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_surface.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix_move.h"
#include "memory.h"
#include "molecule.h"
#include "stl_reader.h"

using namespace LAMMPS_NS;

static constexpr int DELTA = 128;

/* ---------------------------------------------------------------------- */

FixSurface::FixSurface(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

FixSurface::~FixSurface() {}

/* ----------------------------------------------------------------------
   extract lines or tris from a molecule template ID for one or more molecules
   concatenate into single list of points and lines or tris
   identify unique points using hash
------------------------------------------------------------------------- */

void FixSurface::extract_from_molecule(char *molID,
                                       std::map<std::tuple<double,double,double,int>,int> *hash,
                                       int &npoints, int &maxpoints,
                                       Point *&points, int &nlines, Line *&lines,
                                       int &ntris, Tri *&tris)
{
  int dimension = domain->dimension;

  int imol = atom->find_molecule(molID);
  if (imol == -1)
    error->all(FLERR,"Molecule template ID for fix surface does not exist");

  // loop over one or more molecules in molID

  Molecule **onemols = &atom->molecules[imol];
  int nmol = onemols[0]->nset;

  for (int m = 0; m < nmol; m++) {
    if (dimension == 2)
      if (onemols[m]->lineflag == 0)
        error->all(FLERR,"Fix surface molecule must have lines");
    if (dimension == 3)
      if (onemols[m]->triflag == 0)
        error->all(FLERR,"Fix surface molecule must have triangles");

    int nl = onemols[m]->nlines;
    int nt = onemols[m]->ntris;

    nlines += nl;
    ntris += nt;
    lines = (Line *) memory->srealloc(lines,nlines*sizeof(Line),
                                      "surface:lines");
    tris = (Tri *) memory->srealloc(tris,ntris*sizeof(Tri),
                                    "surface:tris");

    // offset line/tri index lists by previous npoints
    // pi,p2,p3 are C-style indices into points vector

    if (dimension == 2) {
      int *molline = onemols[m]->molline;
      int *typeline = onemols[m]->typeline;
      double **epts = onemols[m]->lines;
      int iline = nlines - nl;

      for (int i = 0; i < nl; i++) {
        lines[iline].mol = molline[i];
        lines[iline].type = typeline[i];

        // only lines in the same molecule are connected
        auto key = std::make_tuple(epts[i][0],epts[i][1],0.0,molline[i]);
        if (hash->find(key) == hash->end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface:points");
          }
          (*hash)[key] = npoints;
          points[npoints].x[0] = epts[i][0];
          points[npoints].x[1] = epts[i][1];
          points[npoints].x[2] = 0.0;
          lines[iline].p1 = npoints;
          npoints++;
        } else lines[iline].p1 = (*hash)[key];

        key = std::make_tuple(epts[i][2],epts[i][3],0.0,molline[i]);
        if (hash->find(key) == hash->end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface:points");
          }
          (*hash)[key] = npoints;
          points[npoints].x[0] = epts[i][2];
          points[npoints].x[1] = epts[i][3];
          points[npoints].x[2] = 0.0;
          lines[iline].p2 = npoints;
          npoints++;
        } else lines[iline].p2 = (*hash)[key];

        iline++;
      }
    }

    if (dimension == 3) {
      int *moltri = onemols[m]->moltri;
      int *typetri = onemols[m]->typetri;
      double **cpts = onemols[m]->tris;
      int itri = ntris - nt;

      for (int i = 0; i < nt; i++) {
        tris[itri].mol = moltri[i];
        tris[itri].type = typetri[i];

        // only tris in the same molecule are connected
        auto key = std::make_tuple(cpts[i][0],cpts[i][1],cpts[i][2],moltri[i]);
        if (hash->find(key) == hash->end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface:points");
          }
          (*hash)[key] = npoints;
          points[npoints].x[0] = cpts[i][0];
          points[npoints].x[1] = cpts[i][1];
          points[npoints].x[2] = cpts[i][2];
          tris[itri].p1 = npoints;
          npoints++;
        } else tris[itri].p1 = (*hash)[key];

        key = std::make_tuple(cpts[i][3],cpts[i][4],cpts[i][5],moltri[i]);
        if (hash->find(key) == hash->end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface:points");
          }
          (*hash)[key] = npoints;
          points[npoints].x[0] = cpts[i][3];
          points[npoints].x[1] = cpts[i][4];
          points[npoints].x[2] = cpts[i][5];
          tris[itri].p2 = npoints;
          npoints++;
        } else tris[itri].p2 = (*hash)[key];

        key = std::make_tuple(cpts[i][6],cpts[i][7],cpts[i][8],moltri[i]);
        if (hash->find(key) == hash->end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface:points");
          }
          (*hash)[key] = npoints;
          points[npoints].x[0] = cpts[i][6];
          points[npoints].x[1] = cpts[i][7];
          points[npoints].x[2] = cpts[i][8];
          tris[itri].p3 = npoints;
          npoints++;
        } else tris[itri].p3 = (*hash)[key];

        itri++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   extract triangles from an STL file, can be text or binary
   concatenate into single list of points and tris
   identify unique points using hash
------------------------------------------------------------------------- */

void FixSurface::extract_from_stlfile(char *filename, int stype, int smol,
                                      std::map<std::tuple<double,double,double,int>,int> *hash,
                                      int &npoints, int &maxpoints,
                                      Point *&points, int &ntris, Tri *&tris)
{
  if (domain->dimension == 2)
    error->all(FLERR, "Fix surface cannot use an STL file for 2d simulations");

  if (stype < 1)
    error->all(FLERR, "STL surface type must be >= 1");

  // read tris from STL file
  // stltris = tri coords internal to STL reader

  STLReader *stl = new STLReader(lmp);
  double **stltris;
  int ntris_old = ntris;
  int ntris_new = stl->read_file(filename,stltris);
  ntris += ntris_new;

  tris = (Tri *) memory->srealloc(tris,ntris*sizeof(Tri),"surface:tris");

  // loop over STL tris
  // populate points and tris data structs
  // for each tri: set molID = 1 and type = stype

  for (int itri_new = 0; itri_new < ntris_new; itri_new++) {
    int itri = itri_new + ntris_old;
    tris[itri].mol = smol;
    tris[itri].type = stype;

    // only tris in the same molecule are connected
    auto key = std::make_tuple(stltris[itri_new][0],stltris[itri_new][1],stltris[itri_new][2],smol);
    if (hash->find(key) == hash->end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface:points");
      }
      (*hash)[key] = npoints;
      points[npoints].x[0] = stltris[itri_new][0];
      points[npoints].x[1] = stltris[itri_new][1];
      points[npoints].x[2] = stltris[itri_new][2];
      tris[itri].p1 = npoints;
      npoints++;
    } else tris[itri].p1 = (*hash)[key];

    key = std::make_tuple(stltris[itri_new][3],stltris[itri_new][4],stltris[itri_new][5],smol);
    if (hash->find(key) == hash->end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface:points");
      }
      (*hash)[key] = npoints;
      points[npoints].x[0] = stltris[itri_new][3];
      points[npoints].x[1] = stltris[itri_new][4];
      points[npoints].x[2] = stltris[itri_new][5];
      tris[itri].p2 = npoints;
      npoints++;
    } else tris[itri].p2 = (*hash)[key];

    key = std::make_tuple(stltris[itri_new][6],stltris[itri_new][7],stltris[itri_new][8],smol);
    if (hash->find(key) == hash->end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface:points");
      }
      (*hash)[key] = npoints;
      points[npoints].x[0] = stltris[itri_new][6];
      points[npoints].x[1] = stltris[itri_new][7];
      points[npoints].x[2] = stltris[itri_new][8];
      tris[itri].p3 = npoints;
      npoints++;
    } else tris[itri].p3 = (*hash)[key];
  }

  // delete STL reader

  delete stl;
}

/* ----------------------------------------------------------------------
   create and initialize Connect2d info for global lines
   global lines were read from molecule file(s)
   creates connect2d data structs
------------------------------------------------------------------------- */

void FixSurface::connectivity2d_global(int npoints, int nlines, Line *lines,
                                       Connect2d *&connect2d,
                                       int **&neigh_p1, int **&neigh_p2)
{
  connect2d = (Connect2d *)
    memory->smalloc(nlines*sizeof(Connect2d),"surface:connect2d");

  // setup line end point connectivity lists
  // counts = # of lines containing each end point (including self)
  // plines = ragged 2d array with indices of lines which contain each point

  int *counts;
  memory->create(counts,npoints,"surface:counts");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < nlines; i++) {
    counts[lines[i].p1]++;
    counts[lines[i].p2]++;
  }

  int **plines;
  memory->create_ragged(plines,npoints,counts,"surface:plines");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < nlines; i++) {
    plines[lines[i].p1][counts[lines[i].p1]++] = i;
    plines[lines[i].p2][counts[lines[i].p2]++] = i;
  }

  // p12_counts = # of lines connecting to endpoints p12 of each line
  // do NOT include self

  int *p1_counts,*p2_counts;
  memory->create(p1_counts,nlines,"surface:p1_counts");
  memory->create(p2_counts,nlines,"surface:p2_counts");

  for (int i = 0; i < nlines; i++) {
    p1_counts[i] = counts[lines[i].p1] - 1;
    p2_counts[i] = counts[lines[i].p2] - 1;
  }

  // allocate all ragged arrays which Connect2d will point to

  memory->create_ragged(neigh_p1,nlines,p1_counts,"surface:neigh_p1");
  memory->create_ragged(neigh_p2,nlines,p2_counts,"surface:neigh_p2");

  // set connect2d vector ptrs to rows of corresponding ragged arrays

  for (int i = 0; i < nlines; i++) {
    connect2d[i].np1 = p1_counts[i];
    if (connect2d[i].np1 == 0) connect2d[i].neigh_p1 = nullptr;
    else connect2d[i].neigh_p1 = neigh_p1[i];

    connect2d[i].np2 = p2_counts[i];
    if (connect2d[i].np2 == 0) connect2d[i].neigh_p2 = nullptr;
    else connect2d[i].neigh_p2 = neigh_p2[i];
  }

  // initialize connect2d neigh vectors for each end point of each line
  // do NOT include self

  int j,m;

  for (int i = 0; i < nlines; i++) {
    if (p1_counts[i]) {
      j = 0;
      for (m = 0; m <= p1_counts[i]; m++) {
        if (plines[lines[i].p1][m] == i) continue;
        connect2d[i].neigh_p1[j] = plines[lines[i].p1][m];
        j++;
      }
    }
    if (p2_counts[i]) {
      j = 0;
      for (m = 0; m <= p2_counts[i]; m++) {
        if (plines[lines[i].p2][m] == i) continue;
        connect2d[i].neigh_p2[j] = plines[lines[i].p2][m];
        j++;
      }
    }
  }

  // deallocate counts, plines, p12_counts

  memory->destroy(counts);
  memory->destroy(plines);
  memory->destroy(p1_counts);
  memory->destroy(p2_counts);
}

/* ----------------------------------------------------------------------
   create and initialize Connect3d info for global triangles
   global triangles were read from molecule or STL file(s)
   creates connect3d data structs
------------------------------------------------------------------------- */

int FixSurface::connectivity3d_global(int npoints, int ntris, Tri *tris,
                                      Connect3d *&connect3d,
                                      int **&neigh_e1, int **&neigh_e2,
                                      int **&neigh_e3,
                                      int **&neigh_c1, int **&neigh_c2,
                                      int **&neigh_c3)
{
  int p1,p2,p3;

  connect3d = (Connect3d *)
    memory->smalloc(ntris*sizeof(Connect3d),"surface:connect3d");

  // create hash = map of unique edges
  //   key = <p1,p2> indices of 2 points, in either order
  //   value = index of the unique edge (0 to Nedge-1)
  // tri2edges[i][j] = index of unique edge for tri I and edge J
  // nedges = total count of unique edges, returned

  int **tri2edge;
  memory->create(tri2edge,ntris,3,"surface::tri2edge");

  std::map<std::tuple<int,int>,int> hash;
  int nedges = 0;

  for (int i = 0; i < ntris; i++) {
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;

    auto key1 = std::make_tuple(p1,p2);
    auto key2 = std::make_tuple(p2,p1);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][0] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][0] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][0] = hash[key2];

    key1 = std::make_tuple(p2,p3);
    key2 = std::make_tuple(p3,p2);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][1] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][1] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][1] = hash[key2];

    key1 = std::make_tuple(p3,p1);
    key2 = std::make_tuple(p1,p3);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][2] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][2] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][2] = hash[key2];
  }

  // setup tri edge connectivity lists
  // counts = # of tris containing each edge (including self)
  // etris = ragged 2d array with indices of tris which contain each edge

  int *counts;
  memory->create(counts,nedges,"surface:count");

  for (int i = 0; i < nedges; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    counts[tri2edge[i][0]]++;
    counts[tri2edge[i][1]]++;
    counts[tri2edge[i][2]]++;
  }

  int **etris;
  memory->create_ragged(etris,nedges,counts,"surface:etris");

  for (int i = 0; i < nedges; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    etris[tri2edge[i][0]][counts[tri2edge[i][0]]++] = i;
    etris[tri2edge[i][1]][counts[tri2edge[i][1]]++] = i;
    etris[tri2edge[i][2]][counts[tri2edge[i][2]]++] = i;
  }

  // e123_counts = # of edges connecting to edges e123 of each tri
  // do NOT include self

  int *e1_counts,*e2_counts,*e3_counts;
  memory->create(e1_counts,ntris,"surface:e1_counts");
  memory->create(e2_counts,ntris,"surface:e2_counts");
  memory->create(e3_counts,ntris,"surface:e3_counts");

  for (int i = 0; i < ntris; i++) {
    e1_counts[i] = counts[tri2edge[i][0]] - 1;
    e2_counts[i] = counts[tri2edge[i][1]] - 1;
    e3_counts[i] = counts[tri2edge[i][2]] - 1;
  }

  // allocate all edge ragged arrays which Connect3d will point to

  memory->create_ragged(neigh_e1,ntris,e1_counts,"surface:neigh_e1");
  memory->create_ragged(neigh_e2,ntris,e2_counts,"surface:neigh_e2");
  memory->create_ragged(neigh_e3,ntris,e3_counts,"surface:neigh_e3");

  // set connect3d edge vector ptrs to rows of corresponding edge ragged arrays

  for (int i = 0; i < ntris; i++) {
    connect3d[i].ne1 = e1_counts[i];
    if (connect3d[i].ne1 == 0) connect3d[i].neigh_e1 = nullptr;
    else connect3d[i].neigh_e1 = neigh_e1[i];

    connect3d[i].ne2 = e2_counts[i];
    if (connect3d[i].ne2 == 0) connect3d[i].neigh_e2 = nullptr;
    else connect3d[i].neigh_e2 = neigh_e2[i];

    connect3d[i].ne3 = e3_counts[i];
    if (connect3d[i].ne3 == 0) connect3d[i].neigh_e3 = nullptr;
    else connect3d[i].neigh_e3 = neigh_e3[i];
  }

  // initialize connect3d edge neigh vectors for each edge of each tri
  // do NOT include self

  int j,m;

  for (int i = 0; i < ntris; i++) {
    if (connect3d[i].ne1) {
      j = 0;
      for (m = 0; m < counts[tri2edge[i][0]]; m++) {
        if (etris[tri2edge[i][0]][m] == i) continue;
        connect3d[i].neigh_e1[j] = etris[tri2edge[i][0]][m];
        j++;
      }
    }
    if (connect3d[i].ne2) {
      j = 0;
      for (m = 0; m < counts[tri2edge[i][1]]; m++) {
        if (etris[tri2edge[i][1]][m] == i) continue;
        connect3d[i].neigh_e2[j] = etris[tri2edge[i][1]][m];
        j++;
      }
    }
    if (connect3d[i].ne3) {
      j = 0;
      for (m = 0; m < counts[tri2edge[i][2]]; m++) {
        if (etris[tri2edge[i][2]][m] == i) continue;
        connect3d[i].neigh_e3[j] = etris[tri2edge[i][2]][m];
        j++;
      }
    }
  }

  // deallocate counts, tri2edge, etris, e123_counts

  memory->destroy(counts);
  memory->destroy(tri2edge);
  memory->destroy(etris);
  memory->destroy(e1_counts);
  memory->destroy(e2_counts);
  memory->destroy(e3_counts);

  // setup tri corner point connectivity lists
  // counts = # of tris containing each corner point (including self)
  // ctris = ragged 2d array with indices of tris which contain each point

  memory->create(counts,npoints,"surface:counts");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    counts[tris[i].p1]++;
    counts[tris[i].p2]++;
    counts[tris[i].p3]++;
  }

  int **ctris;
  memory->create_ragged(ctris,npoints,counts,"surface:ctris");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    ctris[tris[i].p1][counts[tris[i].p1]++] = i;
    ctris[tris[i].p2][counts[tris[i].p2]++] = i;
    ctris[tris[i].p3][counts[tris[i].p3]++] = i;
  }

  // c123_counts = # of tris connecting to corner points c123 of each tri
  // do NOT include self or tris which connect to an edge
  // only include tris which only connect at the corner point

  int *c1_counts,*c2_counts,*c3_counts;
  memory->create(c1_counts,ntris,"surface:c1_counts");
  memory->create(c2_counts,ntris,"surface:c2_counts");
  memory->create(c3_counts,ntris,"surface:c3_counts");

  for (int i = 0; i < ntris; i++) {
    c1_counts[i] = counts[tris[i].p1] - 1;
    c1_counts[i] -= connect3d[i].ne3 + connect3d[i].ne1;
    c2_counts[i] = counts[tris[i].p2] - 1;
    c2_counts[i] -= connect3d[i].ne1 + connect3d[i].ne2;
    c3_counts[i] = counts[tris[i].p3] - 1;
    c3_counts[i] -= connect3d[i].ne2 + connect3d[i].ne3;
  }

  // allocate all corner ragged arrays which Connect3d will point to

  memory->create_ragged(neigh_c1,ntris,c1_counts,"surface:neigh_c1");
  memory->create_ragged(neigh_c2,ntris,c2_counts,"surface:neigh_c2");
  memory->create_ragged(neigh_c3,ntris,c3_counts,"surface:neigh_c3");

  // set connect3d corner vector ptrs to rows of corresponding corner ragged arrays

  for (int i = 0; i < ntris; i++) {
    connect3d[i].nc1 = c1_counts[i];
    if (connect3d[i].nc1) connect3d[i].neigh_c1 = neigh_c1[i];
    else connect3d[i].neigh_c1 = nullptr;

    connect3d[i].nc2 = c2_counts[i];
    if (connect3d[i].nc2) connect3d[i].neigh_c2 = neigh_c2[i];
    else connect3d[i].neigh_c2 = nullptr;

    connect3d[i].nc3 = c3_counts[i];
    if (connect3d[i].nc3) connect3d[i].neigh_c3 = neigh_c3[i];
    else connect3d[i].neigh_c3 = nullptr;
  }

  // initialize connect3d corner neigh vectors for each corner of each tri
  // do NOT include self or tris which connect to an edge
  // only include tris which only connect at the corner point

  int n,medge,skipflag;

  for (int i = 0; i < ntris; i++) {
    if (connect3d[i].nc1) {
      j = 0;
      for (m = 0; m < counts[tris[i].p1]; m++) {
        n = ctris[tris[i].p1][m];
        if (n == i) continue;

        skipflag = 0;
        for (medge = 0; medge < connect3d[i].ne3; medge++)
          if (n == connect3d[i].neigh_e3[medge]) skipflag = 1;
        for (medge = 0; medge < connect3d[i].ne1; medge++)
          if (n == connect3d[i].neigh_e1[medge]) skipflag = 1;
        if (skipflag) continue;

        connect3d[i].neigh_c1[j] = n;
        j++;
      }
    }
    if (connect3d[i].nc2) {
      j = 0;
      for (m = 0; m < counts[tris[i].p2]; m++) {
        n = ctris[tris[i].p2][m];
        if (n == i) continue;

        skipflag = 0;
        for (medge = 0; medge < connect3d[i].ne1; medge++)
          if (n == connect3d[i].neigh_e1[medge]) skipflag = 1;
        for (medge = 0; medge < connect3d[i].ne2; medge++)
          if (n == connect3d[i].neigh_e2[medge]) skipflag = 1;
        if (skipflag) continue;

        connect3d[i].neigh_c2[j] = n;
        j++;
      }
    }
    if (connect3d[i].nc3) {
      j = 0;
      for (m = 0; m < counts[tris[i].p3]; m++) {
        n = ctris[tris[i].p3][m];
        if (n == i) continue;

        skipflag = 0;
        for (medge = 0; medge < connect3d[i].ne2; medge++)
          if (n == connect3d[i].neigh_e2[medge]) skipflag = 1;
        for (medge = 0; medge < connect3d[i].ne3; medge++)
          if (n == connect3d[i].neigh_e3[medge]) skipflag = 1;
        if (skipflag) continue;

        connect3d[i].neigh_c3[j] = n;
        j++;
      }
    }
  }

  // deallocate counts, ctris, c123_counts

  memory->destroy(counts);
  memory->destroy(ctris);
  memory->destroy(c1_counts);
  memory->destroy(c2_counts);
  memory->destroy(c3_counts);

  // return edge count

  return nedges;
}
