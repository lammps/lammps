#!/usr/bin/env python

# Script:  refine.py
# Purpose: refine a set of granular surface elements for input to LAMMPS
#          elements in 2d = line segments, elements in 3d = triangles
# Syntax:  refine.py dim insource infile thresh outsource outfile
#          dim = 2 or 3
#          insource = mol for molecule file, stl for STL file
#          infile = molecule or STL filename
#          thresh = max size of resulting surface eleents
#          outsource = mol for molecule file, stl for STL file
#          outfile = molecule or STL filename
# Author:  Steve Plimpton (Sandia), sjplimp at gmail.com

# TODO:
# print stats at each stage
# could also work with local surfs via data and/or dump file
# howto assign tri type, molID for output into molfile
#   maybe should be persisted from input molfile (but not input STL)

import sys
from operator import itemgetter
from math import sqrt

# ----------
# classes
# ----------

class Point:
  def __init__(self,x,y,z):
    self.x = (x,y,z)

class Line:
  def __init__(self,point1,point2):
    self.active = 1
    self.point1 = point1
    self.point2 = point2
    self.p1 = self.p2 = -1

class Edge:                  # unique tri edges
  def __init__(self,p1,p2):
    self.active = 1
    self.p1 = p1
    self.p2 = p2
    self.neigh = []
    self.which = []

class Tri:
  def __init__(self,point1,point2,point3):
    self.active = 1
    self.point1 = point1
    self.point2 = point2
    self.point3 = point3
    self.p1 = self.p2 = self.p3 = -1

# ----------
# functions
# ----------

def error(txt=None):
  if not txt: print("Syntax: refine.py dim insource infile thresh outsource outfile")
  else: print("ERROR:",txt)
  sys.exit()

def distance(x1,x2):
  dx = x1[0] - x2[0]
  dy = x1[1] - x2[1]
  dz = x1[2] - x2[2]
  r = sqrt(dx*dx + dy*dy + dz*dz)
  return r

def midpt(x1,x2):
  xnew = 0.5 * (x1[0] + x2[0])
  ynew = 0.5 * (x1[1] + x2[1])
  znew = 0.5 * (x1[2] + x2[2])
  return (xnew,ynew,znew)

# mybisect in place of Python 3.10 bisect module
# NOTE: explain why Python module does not work
# forward: one > all values before its insertion point
# reverse: one > all values after its insertion point

def mybisect(vec,one,func,reverse):
  value = func(one)
  ilo = 0
  ihi = len(vec)
  niter = 0
  while 1:
    if niter > 10: break
    niter += 1
    imid = (ilo+ihi) // 2
    if ilo == ihi: break
    if not reverse:
      if value <= func(vec[imid]): ihi = imid
      else: ilo = imid+1
    else:
      if value > func(vec[imid]): ihi = imid
      else: ilo = imid+1
  return imid

# debug functions

def print_lines(txt,lines):
  nactive = 0
  for line in lines:
    if line.active: nactive += 1
  print(txt,len(lines),nactive)
  for iline,line in enumerate(lines):
    if not line.active: continue
    print(iline,line.active,distance(points[line.p1].x,points[line.p2].x),
              line.point1,line.point2)

def print_edges(txt,edges):
  nactive = 0
  for edge in edges:
    if edge.active: nactive += 1
  print(txt,len(edges),nactive)
  for iedge,edge in enumerate(edges):
    if not edge.active: continue
    print(iedge,edge.active,distance(points[edge.p1].x,points[edge.p2].x),
          edge.p1,edge.p2,points[edge.p1].x,points[edge.p2].x,
          edge.neigh,edge.which)

def print_tris(txt,tris):
  nactive = 0
  for tri in tris:
    if tri.active: nactive += 1
  print(txt,len(tris),nactive)
  for itri,tri in enumerate(tris):
    if not tri.active: continue
    print(itri,tri.active,distance(points[tri.p1].x,points[tri.p2].x),
              distance(points[tri.p2].x,points[tri.p3].x),
              distance(points[tri.p3].x,points[tri.p1].x),
              tri.point1,tri.point2,tri.point3)

# read a molecule file with lines

def read_molfile_2d(filename):
  with open(filename,"r") as f:
    filelines = f.readlines()[1:]

  nlines = 0

  # read header of molecule file

  for i,oneline in enumerate(filelines):
    oneline = oneline.strip()
    if not oneline: continue
    if oneline[0] == '#': continue
    if nlines: break
    words = oneline.split()
    if words[1] == "lines": nlines = int(words[0])
    if nlines <= 0: error("Molecule file header is invalid for 2d lines")

  # read Lines section of molecule file

  section_word = 0
  count = 0
  lines = []

  for oneline in filelines[i:]:
    oneline = oneline.strip()
    if not oneline: continue
    if section_word == 0:
      if oneline != "Lines":
        error("Molecule file section is invalid for 2d lines")
      section_word = 1
      continue
    count += 1
    words = oneline.split()
    if int(words[0]) != count:
      error("Molecule file section is invalid for 2d lines")
    line = Line((float(words[3]),float(words[4]),0.0),
                (float(words[5]),float(words[6]),0.0))
    lines.append(line)

  if count != nlines:
    error("Molecule file section is invalid for 2d lines")

  return lines

# read a molecule file with tris

def read_molfile_3d(filename):
  with open(filename,"r") as f:
    filelines = f.readlines()[1:]

  ntris = 0

  # read header of molecule file

  for i,oneline in enumerate(filelines):
    oneline = oneline.strip()
    if not oneline: continue
    if oneline[0] == '#': continue
    if ntris: break
    words = oneline.split()
    if words[1] == "triangles": ntris = int(words[0])
    if ntris <= 0: error("Molecule file header is invalid for 3d triangles")

  # read Lines section of molecule file

  section_word = 0
  count = 0
  tris = []

  for oneline in filelines[i:]:
    oneline = oneline.strip()
    if not oneline: continue
    if section_word == 0:
      if oneline != "Triangles":
        error("Molecule file section is invalid for 3d triangles")
      section_word = 1
      continue
    count += 1
    words = oneline.split()
    if int(words[0]) != count:
      error("Molecule file section is invalid for 3d triangles")
    tri = Tri((float(words[3]),float(words[4]),float(words[5])),
              (float(words[6]),float(words[7]),float(words[8])),
              (float(words[9]),float(words[10]),float(words[11])))
    tris.append(tri)

  if count != ntris:
    error("Molecule file section is invalid for 3d triangles")

  return tris

# read an STL file with tris

def read_stlfile(filename):
  with open(filename, "r") as f:
    filelines = f.readlines()

  if not filelines[0].startswith("solid "):
    error("STL file first line is invalid")
  if not filelines[-1].startswith("endsolid "):
    error("STL file last line is invalid")

  filelines = filelines[1:-1]
  if (len(filelines) % 7) != 0:
    error("STL file facet format is invalid")

  tris = []

  iline = 2
  while iline < len(filelines):
    v1 = filelines[iline].split()
    v2 = filelines[iline+1].split()
    v3 = filelines[iline+2].split()
    if v1[0] != "vertex" or v2[0] != "vertex" or v3[0] != "vertex":
      error("STL file vertex format is invalid")
    tri = Tri((float(v1[1]),float(v1[2]),float(v1[3])),
              (float(v2[1]),float(v2[2]),float(v2[3])),
              (float(v3[1]),float(v3[2]),float(v3[3])))
    tris.append(tri)
    iline += 7

  return tris

# create list of unique line end points

def extract_points_2d(lines):
  points = []
  hash = {}
  for line in lines:
    if line.point1 not in hash:
      line.p1 = len(points)
      hash[line.point1] = len(points)
      pt = line.point1
      point = Point(pt[0],pt[1],0.0)
      points.append(point)
    else: line.p1 = hash[line.point1]
    if line.point2 not in hash:
      line.p2 = len(points)
      hash[line.point2] = len(points)
      pt = line.point2
      point = Point(pt[0],pt[1],0.0)
      points.append(point)
    else: line.p2 = hash[line.point2]
  return points

# create list of unique tri corner points

def extract_points_3d(tris):
  points = []
  hash = {}
  for tri in tris:
    if tri.point1 not in hash:
      tri.p1 = len(points)
      hash[tri.point1] = len(points)
      pt = tri.point1
      point = Point(pt[0],pt[1],pt[2])
      points.append(point)
    else: tri.p1 = hash[tri.point1]
    if tri.point2 not in hash:
      tri.p2 = len(points)
      hash[tri.point2] = len(points)
      pt = tri.point2
      point = Point(pt[0],pt[1],pt[2])
      points.append(point)
    else: tri.p2 = hash[tri.point2]
    if tri.point3 not in hash:
      tri.p3 = len(points)
      hash[tri.point3] = len(points)
      pt = tri.point3
      point = Point(pt[0],pt[1],pt[2])
      points.append(point)
    else: tri.p3 = hash[tri.point3]
  return points

# create list of unique tri edges
# also create ehash so can search for edges in refine_3d()
#   key = (p1,p2), value = index into list of edges

def edges_3d(points,tris):
  hash = {}
  for i,tri in enumerate(tris):
    if (tri.p1,tri.p2) in hash: hash[(tri.p1,tri.p2)].append((i,0))
    elif (tri.p2,tri.p1) in hash: hash[(tri.p2,tri.p1)].append((i,0))
    else: hash[(tri.p1,tri.p2)] = [(i,0)]
    if (tri.p2,tri.p3) in hash: hash[(tri.p2,tri.p3)].append((i,1))
    elif (tri.p3,tri.p2) in hash: hash[(tri.p3,tri.p2)].append((i,1))
    else: hash[(tri.p2,tri.p3)] = [(i,1)]
    if (tri.p3,tri.p1) in hash: hash[(tri.p3,tri.p1)].append((i,2))
    elif (tri.p1,tri.p3) in hash: hash[(tri.p1,tri.p3)].append((i,2))
    else: hash[(tri.p3,tri.p1)] = [(i,2)]

  edges = []
  for key,value in hash.items():
    edge = Edge(key[0],key[1])
    edge.neigh = [item[0] for item in value]
    edge.which = [item[1] for item in value]
    edges.append(edge)

  ehash = {}
  for iedge,edge in enumerate(edges):
    ehash[(edge.p1,edge.p2)] = iedge

  return edges,ehash

# successively refine lines which are too long

def refine_2d(points,lines):
  sizes = []
  for iline,line in enumerate(lines):
    sizes.append((iline,distance(points[line.p1].x,points[line.p2].x)))
  sizes = sorted(sizes,key=itemgetter(1),reverse=True)

  # loop until no line is too long

  while sizes[0][1] > thresh:
    #print(sizes)
    iline = sizes[0][0]
    dist = sizes[0][1]

    # add new point

    middle = midpt(points[lines[iline].p1].x,points[lines[iline].p2].x)
    point = Point(middle[0],middle[1],0.0)
    points.append(point)
    npoints = len(points)

    # mark split line as inactive
    # add 2 new lines from split line

    lines[iline].active = 0

    newline = Line(lines[iline].point1,middle)
    newline.p1 = lines[iline].p1
    newline.p2 = npoints - 1
    lines.append(newline)
    inewline1 = len(lines) - 1

    newline = Line(lines[iline].point2,middle)
    newline.p1 = lines[iline].p2
    newline.p2 = npoints - 1
    lines.append(newline)
    inewline2 = len(lines) - 1

    # remove split line from sorted list
    # add 2 new lines to sorted list in appropriate locations

    del sizes[0]

    newsize = (inewline1,0.5*dist)
    isize = mybisect(sizes,newsize,sizeindex,1)
    sizes.insert(isize,newsize)

    newsize = (inewline2,0.5*dist)
    isize = mybisect(sizes,newsize,sizeindex,1)
    sizes.insert(isize,newsize)

# successively refine tris while any edge is too long

def refine_3d(points,edges,ehash,tris):
  sizes = []
  for iedge,edge in enumerate(edges):
    sizes.append((iedge,distance(points[edge.p1].x,points[edge.p2].x)))
  sizes = sorted(sizes,key=itemgetter(1),reverse=True)

  # loop until no tri edge is too long

  while sizes[0][1] > thresh:
    iedge = sizes[0][0]
    edges[iedge].active = 0

    # add new point

    middle = midpt(points[edges[iedge].p1].x,points[edges[iedge].p2].x)
    point = Point(middle[0],middle[1],middle[2])
    points.append(point)
    npoints = len(points)

    # loop over tris sharing the split edge
    #   mark each tri inactive
    #   add 2 new tris to replace inactive tri
    #   add new edge between pair of new tris
    #   reset neighs of non-split and non-new edge for each tri in pair
    # prev_nedges/ntris = number of edges and tris before adding new ones

    neigh = edges[iedge].neigh
    which = edges[iedge].which
    nedges_prev = len(edges)
    ntris_prev = len(tris)

    for itri,iwhich in zip(neigh,which):
      tris[itri].active = 0

      # NOTE:are not updating neighs in edge correctly

      if iwhich == 0:
        newtri = Tri(tris[itri].point1,middle,tris[itri].point3)
        newtri.p1 = tris[itri].p1
        newtri.p2 = npoints - 1
        newtri.p3 = tris[itri].p3
        tris.append(newtri)

        if (tris[itri].p3,tris[itri].p1) in ehash:
          jedge = ehash[(tris[itri].p3,tris[itri].p1)]
        else: jedge = ehash[(tris[itri].p1,tris[itri].p3)]
        index = edges[jedge].neigh.index(itri)
        edges[jedge].neigh[index] = len(tris)-1
        edges[jedge].which[index] = 2

        newtri = Tri(tris[itri].point2,middle,tris[itri].point3)
        newtri.p1 = tris[itri].p2
        newtri.p2 = npoints - 1
        newtri.p3 = tris[itri].p3
        tris.append(newtri)

        if (tris[itri].p2,tris[itri].p3) in ehash:
          jedge = ehash[(tris[itri].p2,tris[itri].p3)]
        else: jedge = ehash[(tris[itri].p3,tris[itri].p2)]
        index = edges[jedge].neigh.index(itri)
        edges[jedge].neigh[index] = len(tris)-1
        edges[jedge].which[index] = 2

        newedge = Edge(tris[itri].p3,npoints-1)
        newedge.neigh = [len(tris)-2,len(tris)-1]
        newedge.which = [1,1]
        edges.append(newedge)
        ehash[(newedge.p1,newedge.p2)] = len(edges)-1

      elif iwhich == 1:
        newtri = Tri(tris[itri].point2,middle,tris[itri].point1)
        newtri.p1 = tris[itri].p2
        newtri.p2 = npoints - 1
        newtri.p3 = tris[itri].p1
        tris.append(newtri)

        if (tris[itri].p1,tris[itri].p2) in ehash:
          jedge = ehash[(tris[itri].p1,tris[itri].p2)]
        else: jedge = ehash[(tris[itri].p2,tris[itri].p1)]
        index = edges[jedge].neigh.index(itri)
        edges[jedge].neigh[index] = len(tris)-1
        edges[jedge].which[index] = 2

        newtri = Tri(tris[itri].point3,middle,tris[itri].point1)
        newtri.p1 = tris[itri].p3
        newtri.p2 = npoints - 1
        newtri.p3 = tris[itri].p1
        tris.append(newtri)

        if (tris[itri].p3,tris[itri].p1) in ehash:
          jedge = ehash[(tris[itri].p3,tris[itri].p1)]
        else: jedge = ehash[(tris[itri].p1,tris[itri].p3)]
        index = edges[jedge].neigh.index(itri)
        edges[jedge].neigh[index] = len(tris)-1
        edges[jedge].which[index] = 2

        newedge = Edge(tris[itri].p1,npoints-1)
        newedge.neigh = [len(tris)-2,len(tris)-1]
        newedge.which = [1,1]
        edges.append(newedge)
        ehash[(newedge.p1,newedge.p2)] = len(edges)-1

      elif iwhich == 2:
        newtri = Tri(tris[itri].point3,middle,tris[itri].point2)
        newtri.p1 = tris[itri].p3
        newtri.p2 = npoints - 1
        newtri.p3 = tris[itri].p2
        tris.append(newtri)

        if (tris[itri].p2,tris[itri].p3) in ehash:
          jedge = ehash[(tris[itri].p2,tris[itri].p3)]
        else: jedge = ehash[(tris[itri].p3,tris[itri].p2)]
        index = edges[jedge].neigh.index(itri)
        edges[jedge].neigh[index] = len(tris)-1
        edges[jedge].which[index] = 2

        newtri = Tri(tris[itri].point1,middle,tris[itri].point2)
        newtri.p1 = tris[itri].p1
        newtri.p2 = npoints - 1
        newtri.p3 = tris[itri].p2
        tris.append(newtri)

        if (tris[itri].p1,tris[itri].p2) in ehash:
          jedge = ehash[(tris[itri].p1,tris[itri].p2)]
        else: jedge = ehash[(tris[itri].p2,tris[itri].p1)]
        index = edges[jedge].neigh.index(itri)
        edges[jedge].neigh[index] = len(tris)-1
        edges[jedge].which[index] = 2

        newedge = Edge(tris[itri].p2,npoints-1)
        newedge.neigh = [len(tris)-2,len(tris)-1]
        newedge.which = [1,1]
        edges.append(newedge)
        ehash[(newedge.p1,newedge.p2)] = len(edges)-1

    # add 2 additional new edges to replace split edge
    # neigh of each new edge depends on order of points in original tri
    # both new edges are first edge 0 in all new tris

    newedge1 = Edge(edges[iedge].p1,npoints-1)
    neigh1 = newedge1.neigh
    newedge2 = Edge(edges[iedge].p2,npoints-1)
    neigh2 = newedge2.neigh

    for itri in range(ntris_prev,len(tris)):
      if edges[iedge].p1 == tris[itri].p1: neigh1.append(itri)
      else: neigh2.append(itri)

    newedge1.which = (len(neigh1))*[0]
    newedge2.which = (len(neigh2))*[0]

    edges.append(newedge1)
    ehash[(newedge1.p1,newedge1.p2)] = len(edges)-1
    edges.append(newedge2)
    ehash[(newedge2.p1,newedge2.p2)] = len(edges)-1

    # remove split edge from sorted list
    # add all new edges to sorted list in appropriate locations

    del sizes[0]

    for iedge in range(nedges_prev,len(edges)):
      newsize = (iedge,distance(points[edges[iedge].p1].x,points[edges[iedge].p2].x))
      isize = mybisect(sizes,newsize,sizeindex,1)
      sizes.insert(isize,newsize)

# function called by bisect() to return entry field that list is sorted on

def sizeindex(entry):
  return entry[-1]

# write a new molecule file with lines

def write_molfile_2d(outfile,lines,infile,thresh):
  nactive = 0
  for line in lines:
    if line.active: nactive += 1

  fp = open(outfile,"w")
  print("LAMMPS molecule file created by refine.py from %s with threshold %g\n" %
        (infile,thresh),file=fp)
  print("%d lines\n" % nactive,file=fp)
  print("Lines\n",file=fp)

  count = 0
  for line in lines:
    if not line.active: continue
    count += 1
    print("%d %d %d %g %g %g %g" % (count,1,1,
                                    line.point1[0],line.point1[1],
                                    line.point2[0],line.point2[1]),file=fp)
  fp.close()

# write a new molecule file with tris

def write_molfile_3d(outfile,tris,infile,thresh):
  nactive = 0
  for tri in tris:
    if tri.active: nactive += 1

  fp = open(outfile,"w")
  print("LAMMPS molecule file created by refine.py from %s with threshold %g\n" %
        (infile,thresh),file=fp)
  print("%d triangles\n" % nactive,file=fp)
  print("Triangles\n",file=fp)

  count = 0
  for tri in tris:
    if not tri.active: continue
    count += 1
    print("%d %d %d %g %g %g %g %g %g %g %g %g" %
          (count,1,1,
           tri.point1[0],tri.point1[1],tri.point1[2],
           tri.point2[0],tri.point2[1],tri.point2[2],
           tri.point3[0],tri.point3[1],tri.point3[2]),
          file=fp)
  fp.close()

# write a new STL file with tris

def write_stlfile(outfile,tris,infile,thresh):
  nactive = 0
  for tri in tris:
    if tri.active: nactive += 1

  fp = open(outfile,"w")
  print("solid STL version of %s with threshold %g" %
        (infile,thresh),file=fp)

  for tri in tris:
    if not tri.active: continue
    print("  facet normal 0.0 0.0 0.0",file=fp)
    print("    outer loop",file=fp)
    print("      vertex %g %g %g" %
          (tri.point1[0],tri.point1[1],tri.point1[2]),file=fp)
    print("      vertex %g %g %g" %
          (tri.point2[0],tri.point2[1],tri.point2[2]),file=fp)
    print("      vertex %g %g %g" %
          (tri.point3[0],tri.point3[1],tri.point3[2]),file=fp)
    print("    endloop",file=fp)
    print("  endfacet",file=fp)

  print("endsolid STL version of %s with threshold %g" %
        (infile,thresh),file=fp)
  fp.close()

# ----------
# main program
# ----------

# parse command-line args

args = sys.argv[1:]
narg = len(args)

if narg != 6: error()

dim = int(args[0])
insource = args[1]
infile = args[2]
thresh = float(args[3])
outsource = args[4]
outfile = args[5]

if dim != 2 and dim != 3: error("Involid dim argument")

# read input file
# lines/tris = list of lines or tris

if insource == "mol" and dim == 2:
  lines = read_molfile_2d(infile)
elif insource == "mol" and dim == 3:
  tris = read_molfile_3d(infile)
elif insource == "stl" and dim == 3:
  tris = read_stlfile(infile)
else:
  error("Invalid insource argument")

# points = list of unique points

if dim == 2: points = extract_points_2d(lines)
if dim == 3: points = extract_points_3d(tris)

# edges = list of unique edges for tris
# ehash = enables search for a p1,p2 edge in edges

if dim == 3: edges,ehash = edges_3d(points,tris)

# perform refinement

if dim == 2: refine_2d(points,lines)
if dim == 3: refine_3d(points,edges,ehash,tris)

# write output file

if outsource == "mol" and dim == 2:
  write_molfile_2d(outfile,lines,infile,thresh)
elif outsource == "mol" and dim == 3:
  write_molfile_3d(outfile,tris,infile,thresh)
elif outsource == "stl" and dim == 3:
  write_stlfile(outfile,tris,infile,thresh)
else:
  error("Invalid outsource argument")
