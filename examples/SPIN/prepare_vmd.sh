#!/bin/bash
# example prepare_vmd.sh /home/jtranch/Documents/lammps/src/dump_VSRSV.lammpstrj
# you will get a return file

echo "vmd script for file $1 is preparing..."

timestamp(){
  date +%s
}

TS=$(timestamp)
FILE=view_${TS}.vmd

cat >${FILE} <<EOF
proc vmd_draw_arrow {mol start end} {
  set middle [vecadd \$start [vecscale 0.9 [vecsub \$end \$start]]]
  graphics \$mol cylinder \$start \$middle radius 0.05
  graphics \$mol cone \$middle \$end radius 0.01 color 3
}

proc vmd_draw_vector {args} {
  set usage {"draw vector {x1 y1 z1} {x2 y2 z2} [scale <s>] [resolution <res>] [radius <r>] [filled <yes/no>]"}
  # defaults
  set scale 2.0
  set res 50
  set radius 0.1
  set filled yes

  if {[llength \$args] < 3} {
    error "wrong # args: should be \$usage"
  }
  set mol    [lindex \$args 0]
  set center [lindex \$args 1]
  set vector [lindex \$args 2]
  if {[llength \$center] != 3 || [llength \$vector] != 3} {
    error "wrong type of args: should be \$usage"
  }

  foreach {flag value} [lrange \$args 3 end] {
    switch -glob \$flag {
      scale  {set scale  \$value}
      res*   {set res    \$value}
      rad*   {set radius \$value}
      fill*  {set filled \$value}
      default {error "unknown option '\$flag': should be \$usage" }
    }
  }

  set vechalf [vecscale [expr \$scale * 0.5] \$vector]
  return [list \\
  [graphics \$mol color yellow]\\
  [graphics \$mol cylinder [vecsub \$center \$vechalf]\\
  [vecadd \$center [vecscale 0.7 \$vechalf]] \\
  radius \$radius resolution \$res filled \$filled] \\
  [graphics \$mol color orange]\\
  [graphics \$mol cone [vecadd \$center [vecscale 0.6 \$vechalf]] \\
  [vecadd \$center \$vechalf] radius [expr \$radius * 2.5] \\
  resolution \$res]]
}

proc vmd_draw_spin {args} {
  global molid
  graphics \$molid delete all
  set frame [molinfo \$molid get frame]
  set natoms [molinfo \$molid get numatoms]
  for {set i 0} {\$i < \$natoms} {incr i} {
    set sel [atomselect top "index \$i"]
    set coords [lindex [\$sel get {x y z}] \$molid]
    set velocities [lindex [\$sel get {vx vy vz}] \$molid]
    draw vector \$coords \$velocities
    set uvx [lindex [\$sel get {vx}] \$molid]
    set uvy [lindex [\$sel get {vy}] \$molid]
    set uvz [lindex [\$sel get {vz}] \$molid]
    \$sel set user [vecadd [vecadd [vecscale \$uvy  \$uvy] [vecscale \$uvz  \$uvz] ] [vecscale \$uvx  \$uvx]]
    \$sel set user \$uvy
    #draw vector \$coords {0.0 uvy 0.0}
  }
  #pbc box -color 3
}

proc enable_trace {} {
  global vmd_frame
  trace variable vmd_frame([molinfo top]) w vmd_draw_spin
}

set molid [mol addfile {$1} type {lammpstrj} autobonds off first 0 last -1 step 1 waitfor all]
scale by 0.5
animate style Loop
enable_trace
EOF
echo "$FILE is ready..."
