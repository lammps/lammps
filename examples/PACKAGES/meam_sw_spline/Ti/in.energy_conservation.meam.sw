# bulk Ti lattice

variable        x index 1
variable        y index 1
variable        z index 1

variable        xx equal 20*$x
variable        yy equal 20*$y
variable        zz equal 20*$z

units           metal
atom_style      atomic

variable       a equal 2.28806
variable       covera equal 1.58111
variable       sqrt3 equal sqrt(3.)
variable       theta equal PI/2.
variable       cos_theta equal round(cos(${theta}))
variable       sin_theta equal round(sin(${theta}))
variable       Dx equal 1.
variable       Dy equal ${covera}
variable       Dz equal sqrt(3.)

lattice        custom ${a} a1 1 0 0 a2 0 ${sqrt3} 0 a3 0 0 ${covera} &
               basis 0.0 0.0 0.0 &
               basis 0.5 0.5 0.0 &
               basis 0.5 0.833333 0.5 &
               basis 0.0 0.333333 0.5 &
               orient x  1             0             0 &
               orient y  0  ${cos_theta}  ${sin_theta} &
               orient z  0 -${sin_theta}  ${cos_theta} &
               spacing ${Dx} ${Dy} ${Dz} &
               origin 0.25 0.25 0.15

variable       Xlo equal -round(9./(${a}*${Dx}))
variable       Xhi equal  round(9./(${a}*${Dx}))
variable       Ylo equal -round(9./(${a}*${Dy}))
variable       Yhi equal  round(9./(${a}*${Dy}))
variable       Zlo equal -round(3./(${a}*${Dz}))
variable       Zhi equal  round(3./(${a}*${Dz}))

region         box prism  ${Xlo} ${Xhi} ${Ylo} ${Yhi} ${Zlo} ${Zhi} 0.0 0.0 0.0

boundary       p p p

create_box      1 box
create_atoms    1 box

pair_style      meam/sw/spline
pair_coeff      * * Ti.meam.sw.spline Ti
mass            * 47.90

velocity        all create 300.0 376847 loop geom

neighbor        1.0 bin
neigh_modify    every 1 delay 5 check yes

fix             1 all nve

thermo 1
thermo_style    custom step vol etotal press pxx pyy pxz
thermo_modify   format 2 %14.8f
thermo_modify   format 3 %14.8f
thermo_modify   format 4 %14.8f
thermo_modify   format 5 %14.8f
thermo_modify   format 6 %14.8f
thermo_modify   format 7 %14.8f


timestep        0.002
thermo          10

run             2000

