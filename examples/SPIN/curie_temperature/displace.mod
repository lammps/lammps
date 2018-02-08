# NOTE: This script should not need to be
# modified. See in.elastic for more info.
#
# Find which reference length to use

if "${dir} == 1" then &
   "variable len0 equal ${lx0}" 
if "${dir} == 2" then &
   "variable len0 equal ${ly0}" 
if "${dir} == 3" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 4" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 5" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 6" then &
   "variable len0 equal ${ly0}" 

# Reset box and simulation parameters

clear
box tilt large
read_restart restart.equil
include potential.mod

# Negative deformation

variable delta equal -${up}*${len0}
variable deltaxy equal -${up}*xy
variable deltaxz equal -${up}*xz
variable deltayz equal -${up}*yz
if "${dir} == 1" then &
   "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box"
if "${dir} == 2" then &
   "change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box"
if "${dir} == 3" then &
   "change_box all z delta 0 ${delta} remap units box"
if "${dir} == 4" then &
   "change_box all yz delta ${delta} remap units box"
if "${dir} == 5" then &
   "change_box all xz delta ${delta} remap units box"
if "${dir} == 6" then &
   "change_box all xy delta ${delta} remap units box"

# Relax atoms positions

minimize ${etol} ${ftol} ${maxiter} ${maxeval}

# Obtain new stress tensor
 
variable tmp equal pxx
variable pxx1 equal ${tmp}
variable tmp equal pyy
variable pyy1 equal ${tmp}
variable tmp equal pzz
variable pzz1 equal ${tmp}
variable tmp equal pxy
variable pxy1 equal ${tmp}
variable tmp equal pxz
variable pxz1 equal ${tmp}
variable tmp equal pyz
variable pyz1 equal ${tmp}

# Compute elastic constant from pressure tensor

variable C1neg equal ${d1}
variable C2neg equal ${d2}
variable C3neg equal ${d3}
variable C4neg equal ${d4}
variable C5neg equal ${d5}
variable C6neg equal ${d6}

# Reset box and simulation parameters

clear
box tilt large
read_restart restart.equil
include potential.mod

# Positive deformation

variable delta equal ${up}*${len0}
variable deltaxy equal ${up}*xy
variable deltaxz equal ${up}*xz
variable deltayz equal ${up}*yz
if "${dir} == 1" then &
   "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box"
if "${dir} == 2" then &
   "change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box"
if "${dir} == 3" then &
   "change_box all z delta 0 ${delta} remap units box"
if "${dir} == 4" then &
   "change_box all yz delta ${delta} remap units box"
if "${dir} == 5" then &
   "change_box all xz delta ${delta} remap units box"
if "${dir} == 6" then &
   "change_box all xy delta ${delta} remap units box"

# Relax atoms positions

minimize ${etol} ${ftol} ${maxiter} ${maxeval}

# Obtain new stress tensor
 
variable tmp equal pe
variable e1 equal ${tmp}
variable tmp equal press
variable p1 equal ${tmp}
variable tmp equal pxx
variable pxx1 equal ${tmp}
variable tmp equal pyy
variable pyy1 equal ${tmp}
variable tmp equal pzz
variable pzz1 equal ${tmp}
variable tmp equal pxy
variable pxy1 equal ${tmp}
variable tmp equal pxz
variable pxz1 equal ${tmp}
variable tmp equal pyz
variable pyz1 equal ${tmp}

# Compute elastic constant from pressure tensor

variable C1pos equal ${d1}
variable C2pos equal ${d2}
variable C3pos equal ${d3}
variable C4pos equal ${d4}
variable C5pos equal ${d5}
variable C6pos equal ${d6}

# Combine positive and negative 

variable C1${dir} equal 0.5*(${C1neg}+${C1pos})
variable C2${dir} equal 0.5*(${C2neg}+${C2pos})
variable C3${dir} equal 0.5*(${C3neg}+${C3pos})
variable C4${dir} equal 0.5*(${C4neg}+${C4pos})
variable C5${dir} equal 0.5*(${C5neg}+${C5pos})
variable C6${dir} equal 0.5*(${C6neg}+${C6pos})

# Delete dir to make sure it is not reused

variable dir delete
