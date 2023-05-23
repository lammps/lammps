# Definition of SNAP+ZBL potential.
variable zblcutinner equal 0.1 
variable zblcutouter equal 0.2
variable zblz equal  10 

# Specify hybrid with SNAP, ZBL, and long-range C_SNAP_2021.10.15.quadratic.ulomb

pair_style hybrid/overlay &
zbl ${zblcutinner} ${zblcutouter} &
snap
pair_coeff 1 1 zbl ${zblz} ${zblz}
pair_coeff * * snap C_SNAP_2021.10.15.quadratic.snapcoeff C_SNAP_2021.10.15.quadratic.snapparam C
