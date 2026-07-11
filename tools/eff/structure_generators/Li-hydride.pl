#!/usr/bin/perl

# Usage: ./Li-hydride.pl <nx> <ny> <nz> > data.Li-hydride
#
# Generates a lithium hydride solid.

$nx = shift(@ARGV);
$ny = shift(@ARGV);
$nz = shift(@ARGV);

$L = 6.785;  # eFF optimized (7.72 expt)

# This part changes for different lattices
@xunit = (0, 0.5, 0, 0.5, 0.5, 0, 0, 0.5);
@yunit = (0, 0.5, 0.5, 0, 0, 0.5, 0, 0.5);
@zunit = (0, 0, 0.5, 0.5, 0, 0, 0.5, 0.5);

$r_elec = 0.7;
$r2_elec = 2.2;

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
        $xnuc[$idx] = $x * $Lx + $xunit[$i] * $L + 0.5 * $L;
        $ynuc[$idx] = $y * $Ly + $yunit[$i] * $L + 0.5 * $L;
        $znuc[$idx] = $z * $Lz + $zunit[$i] * $L + 0.5 * $L;
        $idx++;
      } 
    }
  }
} 

$numnuc = $idx;

# Print length of supercell

printf("Created with Li-hydride.pl\n\n");
printf("%d atoms\n",$numnuc+$numnuc+$numnuc);
printf("3 atom types\n\n");
printf("%f %f xlo xhi\n", 0, $Lx * $nx);
printf("%f %f ylo yhi\n", 0, $Ly * $ny);
printf("%f %f zlo zhi\n\n", 0, $Lz * $nz);
printf("Masses\n\n");
printf("1 6.941000\n");
printf("2 1.007940\n");
printf("3 1.000000\n\n");
printf("Atoms\n\n");

$j = 0;
# Print out the nuclei
for ($i = 0; $i < $numnuc; $i += 8)
{
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 1, 3.0, 0, 0.0,$xnuc[$i], $ynuc[$i], $znuc[$i]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 1, 3.0, 0, 0.0,$xnuc[$i+1], $ynuc[$i+1], $znuc[$i+1]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 1, 3.0, 0, 0.0,$xnuc[$i+2], $ynuc[$i+2], $znuc[$i+2]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 1, 3.0, 0, 0.0,$xnuc[$i+3], $ynuc[$i+3], $znuc[$i+3]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 1.0, 0, 0.0,$xnuc[$i+4], $ynuc[$i+4], $znuc[$i+4]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 1.0, 0, 0.0,$xnuc[$i+5], $ynuc[$i+5], $znuc[$i+5]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 1.0, 0, 0.0,$xnuc[$i+6], $ynuc[$i+6], $znuc[$i+6]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 1.0, 0, 0.0,$xnuc[$i+7], $ynuc[$i+7], $znuc[$i+7]);
}

# Print out the electrons
for ($i = 0; $i < $numnuc; $i += 8)
{
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, 1, $r_elec,$xnuc[$i  ], $ynuc[$i  ], $znuc[$i  ]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, -1, $r_elec,$xnuc[$i  ], $ynuc[$i  ], $znuc[$i  ]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, 1, $r_elec,$xnuc[$i+1], $ynuc[$i+1], $znuc[$i+1]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, -1, $r_elec,$xnuc[$i+1], $ynuc[$i+1], $znuc[$i+1]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, 1, $r_elec,$xnuc[$i+2], $ynuc[$i+2], $znuc[$i+2]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, -1, $r_elec,$xnuc[$i+2], $ynuc[$i+2], $znuc[$i+2]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, 1, $r_elec,$xnuc[$i+3], $ynuc[$i+3], $znuc[$i+3]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, -1, $r_elec,$xnuc[$i+3], $ynuc[$i+3], $znuc[$i+3]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, 1, $r_elec,$xnuc[$i+4], $ynuc[$i+4], $znuc[$i+4]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, -1, $r_elec,$xnuc[$i+4], $ynuc[$i+4], $znuc[$i+4]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, 1, $r_elec,$xnuc[$i+5], $ynuc[$i+5], $znuc[$i+5]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, -1, $r_elec,$xnuc[$i+5], $ynuc[$i+5], $znuc[$i+5]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, 1, $r_elec,$xnuc[$i+6], $ynuc[$i+6], $znuc[$i+6]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, -1, $r_elec,$xnuc[$i+6], $ynuc[$i+6], $znuc[$i+6]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, 1, $r_elec,$xnuc[$i+7], $ynuc[$i+7], $znuc[$i+7]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 3, 0.0, 1, $r_elec,$xnuc[$i+7], $ynuc[$i+7], $znuc[$i+7]);
}

