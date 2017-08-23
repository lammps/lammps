#!/usr/bin/env python

err_msg = """
Usage:

   generate_system_lt.py n < monomer_coords.raw > system.lt

Example:

   generate_system_lt.py 30118 47 < coords.raw > system.lt

Explanation:
     ARGUMENTS:
         n = total length of the polymer (in monomers)
         L = the (average) length of each condensin interval (Poisson-
             distributed)  This is also 1/probability that each monomer
             is a "condensin monomer".

   (Note: 30117 ~= 128000/4.25, but using 30118 makes interpolation cleaner,
          and 47 = 200/4.25. Note that 128000 and 200 are for the 10nm model.
          See the supplemental section of Naumova et al Science 2013, p 18.)

"""


import sys
import random
from math import *

# Parse the argument list:
if len(sys.argv) <= 2:
    sys.stderr.write("Error:\n\nTypical Usage:\n\n"+err_msg+"\n")
    exit(1)
N=int(sys.argv[1])
L=float(sys.argv[2])
if len(sys.argv) > 3:
    delta_x = float(sys.argv[3])
else:
    delta_x = 2.0
if len(sys.argv) > 4:
    x_offset = float(sys.argv[4])
else:
    x_offset = -((N-1.0)/2) * delta_x


coords = [[0.0, 0.0, 0.0] for i in range(0,N)]
lines = sys.stdin.readlines()
if len(lines) != N:
    sys.stderr.write("Error: Number of lines in input file ("+str(len(lines))+")\n"
                     "       does not match first argument ("+str(N)+")\n")
    exit(1)
for i in range(0, N):
    coords[i] = list(map(float, lines[i].split()))

# Now calculate the box_boundaries:
box_bounds_min = [0.0, 0.0, 0.0]
box_bounds_max = [0.0, 0.0, 0.0]
for i in range(0, N):
    for d in range(0, 3):
        if i == 0:
            box_bounds_min[d] = coords[i][d]
            box_bounds_max[d] = coords[i][d]
        else:
            if coords[i][d] > box_bounds_max[d]:
                box_bounds_max[d] = coords[i][d]
            if coords[i][d] < box_bounds_min[d]:
                box_bounds_min[d] = coords[i][d]

# Now scale the box boundaries outward by 50%
box_scale = 1.5
for d in range(0,3):
    box_bounds_cen    = 0.5*(box_bounds_max[d] + box_bounds_min[d])
    box_bounds_width  = box_bounds_max[d] - box_bounds_min[d]
    box_bounds_min[d] = box_bounds_cen - 0.5*box_scale*box_bounds_width
    box_bounds_max[d] = box_bounds_cen + 0.5*box_scale*box_bounds_width

# Now calculate the direction each molecule should be pointing at:
direction_vects = [[0.0, 0.0, 0.0] for i in range(0,N)]
for d in range(0, 3):
     direction_vects[0][d] = coords[1][d] - coords[0][d]
     direction_vects[N-1][d] = coords[N-1][d] - coords[N-2][d]
for i in range(1, N-1):
    for d in range(0, 3):
         direction_vects[i][d] = coords[i+1][d] - coords[i-1][d]

# Optional: normalize the direction vectors
for i in range(1, N-1):
    direction_len = 0.0
    for d in range(0, 3):
        direction_len += (direction_vects[i][d])**2
    direction_len = sqrt(direction_len)
    for d in range(0, 3):
        direction_vects[i][d] /= direction_len

# Now, begin writing the text for the system.lt file:

sys.stdout.write(
"""
import "monomer.lt"      # <-- defines "Monomer"
import "condensin.lt"    # <-- defines "CondensinMonomer"


"""
)



# Figure out which monomers are "Monomers" and which monomers are
# "CondensinMonomers"

ic = 0 # count the number of condensins added so far
condensin_is_here = [False for i in range(0, N)]
for i in range(0, N):
    #add_condensin_here = random.random() <  (1.0 / L)
    add_condensin_here = random.random() <  (1.0 / (L-2.0))

    # We do not allow condensin at successive sites separated by less than 2
    # subunits (the "L-2.0" above is to compensate for this)
    if (((i > 0) and condensin_is_here[i-1]) or
        ((i > 1) and condensin_is_here[i-2])):
        add_condensin_here = False

    if add_condensin_here:
        condensin_is_here[i] = True
        ic += 1
Nc = ic


ic = 0
for i in range(0, N):
    if condensin_is_here[i]:
        sys.stdout.write("condensins["+str(ic)+"] = new CondensinMonomer.scale(0.5,0.8,0.8).rotvv(1,0,0,")
        ic+=1
    else:
        sys.stdout.write("monomers["+str(i)+"] = new Monomer.scale(0.5,0.8,0.8).rotvv(1,0,0,")
    sys.stdout.write(str(direction_vects[i][0])+","
                     +str(direction_vects[i][1])+","
                     +str(direction_vects[i][2])+
                     ").move("
                     +str(coords[i][0])+","
                     +str(coords[i][1])+","
                     +str(coords[i][2])+")\n")

    #if condensin_is_here[i]:
    #    if i < N-1:
    #        sys.stdout.write("\n"
    #                         "#(override the dihedral angle for this monomer)\n"
    #                         "write(\"Data Dihedrals\") {\n"
    #                         "  $dihedral:twistor"+str(i+1)+" @dihedral:CondensinMonomer/TWISTOR $atom:monomers["+str(i)+"]/t $atom:monomers["+str(i)+"]/c $atom:monomers["+str(i+1)+"]/c $atom:monomers["+str(i+1)+"]/t\n"
    #                         "}\n"
    #                         "\n")



sys.stdout.write(
"""

# ---------------- simulation box -----------------

# Now define a box big enough to hold a polymer with this (initial) shape

"""
)


sys.stdout.write("write_once(\"Data Boundary\") {\n"
                 +str(box_bounds_min[0])+"  "+str(box_bounds_max[0])+" xlo xhi\n"
                 +str(box_bounds_min[1])+"  "+str(box_bounds_max[1])+" ylo yhi\n"
                 +str(box_bounds_min[2])+"  "+str(box_bounds_max[2])+" zlo zhi\n"
                 "}\n\n\n")


sys.stdout.write(
"""
# What kind of boundary conditions are we using?

write_once("In Init") {
  boundary s s s      # <-- boundary conditions in x y z directions
  #boundary p p p      # <-- boundary conditions in x y z directions
}
# "p" stands for "periodic"
# "s" stands for "shrink-wrapped" (non-periodic)


# ---- Bonds ----


write_once("In Settings") {
  #  10nm model:
  #bond_coeff @bond:backbone harmonic 100.0 1.0
  #  30nm fiber (4.25^(1/3)=1.6198059006387417)
  bond_coeff @bond:backbone harmonic 100.0 1.6198059006387417
}


"""
)


sys.stdout.write("write(\"Data Bonds\") {\n")

# Old bond-loop was simple:
#for i in range(0, N-1):
#     sys.stdout.write("  $bond:b"+str(i+1)+" @bond:backbone $atom:monomers["+str(i)+"]/a  $atom:monomers["+str(i+1)+"]/a\n")

ic = 0
for i in range(0, N-1):
    #sys.stderr.write("i="+str(i)+", ic="+str(ic)+", Nc="+str(Nc)+"\n")

    # Figure out if the first atom in the bond pair
    # belongs to a regular Monomer or a CondensinMonomer
    if condensin_is_here[i]:
        sys.stdout.write("  $bond:b"+str(i+1)+" @bond:backbone $atom:condensins["+str(ic)+"]/a")
        ic+=1
    else:
        sys.stdout.write("  $bond:b"+str(i+1)+" @bond:backbone $atom:monomers["+str(i)+"]/a")

    # Do the same thing for the second atom in the bond pair
    if condensin_is_here[i+1]:
        assert(ic<Nc)
        sys.stdout.write("  $atom:condensins["+str(ic)+"]/a\n")
    else:
        sys.stdout.write("  $atom:monomers["+str(i+1)+"]/a\n")

sys.stdout.write("}\n\n\n")


sys.stdout.write("""

write_once("Data Angles By Type") {
  @angle:backbone  @atom:*  @atom:*  @atom:*  @bond:backbone  @bond:backbone
}

write_once("In Settings") {
  # Most parameters here were taken from the supplemental material of
  # Naumova et al. Science 2013 (simulations by Maxim Imakaev, see Supp Mat)
  #angle_coeff @angle:backbone cosine 5.0                 #<-10nm fiber
  angle_coeff @angle:backbone cosine 1.1764705882352942  #<-30nm fiber
}

""")


sys.stdout.write(
"""

# ---- Condensins randomly located on the polymer ----

# Stage 1:
# Add bonds between consecutive condensin anchors.
# Imakaev calls this "stage 1: linear compaction":

"""
)

sys.stdout.write("write(\"Data Bonds\") {\n")
ic = 0
for i in range(0, N):
    if condensin_is_here[i]:
        if (0 < ic):
            sys.stdout.write("  $bond:bstage1_"+str(ic-1)+"  @bond:stage1  $atom:condensins["+str(ic-1)+"]/a $atom:condensins["+str(ic)+"]/a\n")
        ic += 1

sys.stdout.write("}\n\n")



sys.stdout.write("""

# Stage 2:
# Add additional bonds between all pairs of condensin anchors
# in a window of |ic-jc| <= 30 anchors (along the chain).
# In the paper, they call this stage 2 axial compression.

write("Data Bonds") {
""")

jcwindowsize = 30
for ic in range(0, Nc):
    #jcmin = max(ic-jcwindowsize, 0)
    #jcmax = min(ic+jcwindowsize, Nc-1)
    jcmax = min(ic+jcwindowsize, Nc-1)
    for jc in range(ic+2, jcmax+1):
        sys.stdout.write("  $bond:bstage2_"+str(ic)+"_"+str(jc)+"  @bond:stage2  $atom:condensins["+str(ic)+"]/a  $atom:condensins["+str(jc)+"]/a\n")

sys.stdout.write("}\n")


sys.stdout.write("""

write_once("In Settings") {
  # stage 1 bonds are initially off
  bond_coeff @bond:stage1 harmonic 0.0 0.5   # <--(we can override this later)"
  # stage 2 bonds are initially off
  bond_coeff @bond:stage2 harmonic 0.0 0.0   # <--(we can override this later)"
}

""")

sys.stdout.write("\n\n# "+str(Nc)+" condensin molecules added\n\n")
sys.stderr.write("\n\n# "+str(Nc)+" condensin molecules added\n\n")
