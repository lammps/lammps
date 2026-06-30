units		lj
atom_style	atomic
atom_modify 	map yes sort 0 0

boundary	f p p
lattice		fcc ${rho} origin ${eps} ${eps} ${eps}

region		box block 0 ${nx} 0 ${ny} 0 ${nz} units lattice

create_box	1 box
create_atoms	1 region box
mass		1 ${mass}

# define property to store timestep it crosses threshold velocity

fix             myidump all property/atom i_tshock d_xshock &
d_vxhist010 &
d_vxhist020 &
d_vxhist030 &
d_vxhist040 &
d_vxhist050 &
d_vxhist060 &
d_vxhist070 &
d_vxhist080 &
d_vxhist090 &
d_vxhist100 &
d_vxhist110 &
d_vxhist120 &
d_vxhist130 &
d_vxhist140 &
d_vxhist150 &
d_vxhist160 &
d_vxhist170 &
d_vxhist180 &
d_vxhist190 &
d_vxhist200

set             group all i_tshock -1
set             group all d_xshock 0.0
set             group all &
d_vxhist010 -1e20 &
d_vxhist020 -1e20 &
d_vxhist030 -1e20 &
d_vxhist040 -1e20 &
d_vxhist050 -1e20 &
d_vxhist060 -1e20 &
d_vxhist070 -1e20 &
d_vxhist080 -1e20 &
d_vxhist090 -1e20 &
d_vxhist100 -1e20 &
d_vxhist110 -1e20 &
d_vxhist120 -1e20 &
d_vxhist130 -1e20 &
d_vxhist140 -1e20 &
d_vxhist150 -1e20 &
d_vxhist160 -1e20 &
d_vxhist170 -1e20 &
d_vxhist180 -1e20 &
d_vxhist190 -1e20 &
d_vxhist200 -1e20

velocity	all create ${temp} 87287

pair_style	lj/cubic 
pair_coeff 	* * 1.0 1.0

neighbor	0.3 bin
neigh_modify	every 20 delay 0 check no

thermo_style  	custom step temp pe press lx ly lz 
run 0

timestep	${dt}
variable        tdamp equal ${dt}*10.0
fix		mynvt all nvt temp ${temp} ${temp} ${tdamp} 
fix 		xwalls all wall/reflect xlo EDGE xhi EDGE

run		${nequil}

# nemd simulation

reset_timestep	0

# estimate location of shock front assuming ideal kinematics
# i.e. P(t) = P0 * (L0+u0*t-xs)/L0, requires u1 = 0
compute	   	mymom all momentum
variable	xshock equal (1.0+c_mymom[1]/(atoms*${mass})/${up})*lx-${up}*time
fix             xshockprev all vector ${nthermo} v_xshock nmax 2
variable        ushock equal (f_xshockprev[2]-f_xshockprev[1])/(${nthermo}*${dt})+${up}
run             ${nthermo} # needed to fil both entries in xshockprev vector

thermo          ${nthermo}
thermo_style    custom step temp pe density etotal press v_xshock v_ushock

# add impact velocity in minus x

variable	upneg equal -${up}
velocity	all set ${upneg} 0.0 0.0 sum yes units box

unfix		mynvt
fix		mynve all nve

# create group for bulk atoms, to avoid surface artefacts

variable        bulklo equal xlo+${leftbulkpad}
variable        bulkhi equal xhi-${rightbulkpad}
region		bulk block ${bulklo} ${bulkhi} INF INF INF INF units box
group           bulk region bulk

