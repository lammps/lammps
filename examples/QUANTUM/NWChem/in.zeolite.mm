# MM for SiO2 zeolite with one methane molecule

# CHIK potential
# EPL, Carre, Horbach, Ispas, Kob, 82, 17001 (2008)
# B = 1/rho

#q Si = 1.910418
#q O = -0.955209
#A OO = 659.595398 eV
#B OO = 2.590066 1/Ang
#C OO = 26.836679 eV-Ang^6
#A SiO = 27029.419922 eV
#B SiO = 5.158606 1/Ang
#C SiO = 148.099091 eV-Ang^6
#A SiSi = 3150.462646 eV
#B SiSi = 2.851451 1/Ang
#C SiSi = 626.7519553 eV-Ang^6

# LJ params for methane and O from Table 1
# Bhatia and Nicholson, J Phys Chem C, 2012, 116, 2344-2355.

#q C = -0.66
#Q H = 0.165
#sigma C = 0.34 nm
#sigma H = 0.265 nm
#sigma O = 0.28 nm
#eps/kB C = 55.082 K = 0.004745993 eV
#eps/kB H = 7.905 K = 0.000681113 eV
#eps/kB O = 492.7 K = 0.0424522 eV

# LJ params for silicon
#e-Journal of Surf Sci and Nanotech, Inui and Iwasaki, 15, 40-49 (2017)

#sigma Si = 3.826 Ang
#eps Si = 17.4 meV = 0.0174 eV

# C-H bond and methane angle params

#OPLS C-H bond k = 29.40 ev/Ang^2
#C-H bond r0 = 1.09 Angs
#methane angles = 109.5 degrees
#C-H angle k/kB = 2000 K/rad^2

# conversions

#1 eV = 11606 K
#1 eV = 23.0609 kcal/mole
#1 kcal/mole = 503.2761 K
#1 kcal = 4.814 kJoule

# -------------------------

units		metal
atom_style	full

bond_style      harmonic
angle_style     harmonic

read_data       data.zeolite

group           mm type 1 2
group           qm type 3 4

# pair style must define stand-alone short-range Coulombics
# arithmetic mixing

pair_style      hybrid/overlay buck 6.5 lj/cut 6.5 coul/cut 6.5

pair_coeff      1 1 buck 3150.462646 0.35032282 626.7519553
pair_coeff      2 2 buck 659.595398 0.38609055 26.836679
pair_coeff      1 2 buck 27029.419922 0.19385082 148.099091
pair_coeff      1 2 buck 27029.419922 0.19385082 148.099091
pair_coeff      1 3 lj/cut 0.009087 3.613
pair_coeff      1 4 lj/cut 0.00344258 3.238
pair_coeff      2 3 lj/cut 0.01419429 3.1
pair_coeff      2 4 lj/cut 0.00537724 2.725
pair_coeff      3 3 lj/cut 0.004746 3.4
pair_coeff      4 4 lj/cut 0.00068111 2.65
pair_coeff      * * coul/cut

bond_style      harmonic
bond_coeff      1 29.40 1.09

angle_style      harmonic
angle_coeff      1 0.172325 109.5

#velocity        all create 300.0 458732

neighbor	1.0 bin
neigh_modify	delay 0 every 1 check yes

# dynamic or frozen zeolite

fix		1 all nve
#fix		1 qm nve

timestep        0.0001

thermo_style    custom step cpu temp ke evdwl ecoul epair emol elong &
                pe etotal press

thermo          1

run		3
