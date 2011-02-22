#!/usr/bin/perl

# Usage: Li-solid.pl <nx> <ny> <nz> > data.Li
#
# Generates a lithium solid LAMMPS data input file
# with FCC nuclear positions and interstitial electron positions.

$nx = shift(@ARGV);
$ny = shift(@ARGV);
$nz = shift(@ARGV);

$L = 4.419688383648;  # eFF optimized (8.32 expt)

# This part changes for different lattices
@xunit = (0, 0.2645886245, 0, 0.2645886245, 0.2645886245, 0, 0, 0.2645886245);
@yunit = (0, 0.2645886245, 0.2645886245, 0, 0, 0.2645886245, 0, 0.2645886245);
@zunit = (0, 0, 0.2645886245, 0.2645886245, 0, 0, 0.2645886245, 0.2645886245);

$r_elec = 0.3704240743;
$r_elec2 = 1.587531747;

$Lx = $L;
$Ly = $L;
$Lz = $L;

$idx = 0;
for ($x = 0; $x < $nx; $x++)
{
 for ($y = 0; $y < $ny; $y++)
 {
   for ($z = 0; $z < $nz; $z++)
   {
     for ($i = 0; $i <= $#xunit; $i++)
     {
       $xnuc[$idx] = $x * $Lx + $xunit[$i] * $L + 0.2645886245 * $L;
       $ynuc[$idx] = $y * $Ly + $yunit[$i] * $L + 0.2645886245 * $L;
       $znuc[$idx] = $z * $Lz + $zunit[$i] * $L + 0.2645886245 * $L;
       $idx++;
     }
   }
 }
} 

$numnuc = $idx;

# Print length of supercell

printf("Created by AJB\n\n");
printf("%d atoms\n",2*$numnuc);
printf("2 atom types\n\n");
printf("%f %f xlo xhi\n", 0, $Lx * $nx);
printf("%f %f ylo yhi\n", 0, $Ly * $ny);
printf("%f %f zlo zhi\n\n", 0, $Lz * $nz);
printf("Masses\n\n");
printf("1 6.941000\n");
printf("2 1.000000\n\n");
printf("Atoms\n\n");

$j = 0;
# Print out the nuclei
for ($i = 0; $i < $numnuc; $i += 8)
{
 printf("%i %i %f %i %f %f %f %f\n", $j+=1,1,3.0, 0, 0.0,$xnuc[$i], $ynuc[$i], $znuc[$i]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1,1,3.0, 0, 0.0,$xnuc[$i+1], $ynuc[$i+1], $znuc[$i+1]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1,1,3.0, 0, 0.0,$xnuc[$i+2], $ynuc[$i+2], $znuc[$i+2]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1,1,3.0, 0, 0.0,$xnuc[$i+3], $ynuc[$i+3], $znuc[$i+3]);
}

# Print out the core electrons
for ($i = 0; $i < $numnuc; $i += 8)
{
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, 1, $r_elec,$xnuc[$i  ], $ynuc[$i  ], $znuc[$i  ]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, -1, $r_elec,$xnuc[$i  ], $ynuc[$i  ], $znuc[$i  ]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, 1, $r_elec,$xnuc[$i+1], $ynuc[$i+1], $znuc[$i+1]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, -1, $r_elec,$xnuc[$i+1], $ynuc[$i+1], $znuc[$i+1]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, 1, $r_elec,$xnuc[$i+2], $ynuc[$i+2], $znuc[$i+2]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, -1, $r_elec,$xnuc[$i+2], $ynuc[$i+2], $znuc[$i+2]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, 1, $r_elec,$xnuc[$i+3], $ynuc[$i+3], $znuc[$i+3]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, -1, $r_elec,$xnuc[$i+3], $ynuc[$i+3], $znuc[$i+3]);
}

# Print out the valence electrons
for ($i = 0; $i < $numnuc; $i += 8)
{
 if (rand() < .5) {$spin = 1;} else {$spin = -1;}  
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, $spin, $r_elec2,$xnuc[$i+4], $ynuc[$i+4], $znuc[$i+4]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, -$spin, $r_elec2,$xnuc[$i+5], $ynuc[$i+5], $znuc[$i+5]);
 if (rand() < .5) {$spin = 1;} else {$spin = -1;}  
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, $spin, $r_elec2,$xnuc[$i+6], $ynuc[$i+6], $znuc[$i+6]);
 printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, -$spin, $r_elec2,$xnuc[$i+7], $ynuc[$i+7], $znuc[$i+7]);
}

