# Ar in metal units

# simulation params in reduced units
# settable from command line
# epsilon, sigma, mass set below

variable	x index 5
variable	y index 5
variable	z index 5
variable        rhostar index 0.8842
variable        dt index 0.005
variable        cutoff index 2.5
variable        skin index 0.3
variable        tinitial index 1.0
variable        nthermo index 10
variable        nsteps index 100

# physical constants from update.cpp

variable        kb index 8.617343e-5          # kB in eV/K
variable        avogadro index 6.02214129e23  # Avogadro's number

# Ar properties in metal units

variable        epskb index 117.7             # LJ epsilon/kB in degrees K
variable        sigma index 3.504             # LJ sigma in Angstroms
variable        epsilon equal ${epskb}*${kb}  # LJ epsilon in eV
variable        mass index 39.95              # mass in g/mole

# scale factors

# sigma = scale factor on distance, converts reduced distance to Angs
# epsilon = scale factor on energy, converts reduced energy to eV
# tmpscale = scale factor on temperature, converts reduced temp to degrees K
# tscale = scale factor on time, converts reduced time to ps
#   formula is t = t* / sqrt(epsilon/mass/sigma^2), but need t in fs
#   use epsilon (Joule), mass (kg/atom), sigma (meter) to get t in seconds
# pscale = scale factor on pressure, converts reduced pressure to bars
#   formula is P = P* / (sigma^3/epsilon), but need P in atmospheres
#   use sigma (meter), epsilon (Joule) to get P in nt/meter^2, convert to bars

variable        eVtoJoule index 1.602e-19     # convert eV to Joules
variable        NtMtoAtm equal 1.0e-5         # convert Nt/meter^2 to bars

variable        tmpscale equal ${epskb}
variable        epsilonJ equal ${epsilon}*${eVtoJoule}
variable        massKgAtom equal ${mass}/1000.0/${avogadro}
variable        sigmaM equal ${sigma}/1.0e10
variable        sigmaMsq equal ${sigmaM}*${sigmaM}
variable        tscale equal 1.0e12/sqrt(${epsilonJ}/${massKgAtom}/${sigmaMsq})
variable        sigmaM3 equal ${sigmaM}*${sigmaM}*${sigmaM}
variable        pscale equal ${NtMtoAtm}/(${sigmaM3}/(${epsilonJ}))

# variables
# alat = lattice constant in Angs (at reduced density rhostar)
# temp = reduced temperature for output
# epair,emol,etotal = reduced epair,emol,etotal energies for output
# press = reduced pressure for output

variable        alat equal (4.0*${sigma}*${sigma}*${sigma}/${rhostar})^(1.0/3.0)
variable        temp equal temp/${tmpscale}
variable        epair equal epair/${epsilon}
variable        emol equal emol/${epsilon}
variable        etotal equal etotal/${epsilon}
variable        press equal press/${pscale}

# same script as in.ar.lj

units		metal
atom_style	atomic

lattice		fcc ${alat}
region		box block 0 $x 0 $y 0 $z
create_box	1 box
create_atoms	1 box
mass		1 ${mass}

velocity	all create $(v_tinitial*v_epskb) 12345

pair_style	lj/cut $(v_cutoff*v_sigma)
pair_coeff	1 1 ${epsilon} ${sigma}

neighbor	$(v_skin*v_sigma) bin
neigh_modify	delay 0 every 20 check no

fix		1 all nve

timestep	$(v_dt*v_tscale)

# columns 2,3,4 = temp,pe,press in metal units
# columns 5-9 = temp,energy.press in reduced units, compare to in.ar.lj
# need to include metal unit output to enable use of reduced variables

thermo_style    custom step temp pe press v_temp v_epair v_emol v_etotal v_press
thermo_modify	norm yes
thermo		${nthermo}

run		${nsteps}
