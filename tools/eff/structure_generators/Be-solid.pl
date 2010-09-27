#!/usr/bin/perl

# Usage: Be-solid.pl <nx> <ny> <nz> > data.Be
#
# Generates a beryllium solid LAMMPS data file.

$nx = shift(@ARGV);
$ny = shift(@ARGV);
$nz = shift(@ARGV);

$La = 4.592;  # eFF optimized (4.33 expt)
$Lc = 7.025;  # eFF optimized (6.78 expt)

# This part changes for different lattices
@xunit = (0, 0.5, 0, 0.5); 
@yunit = (0, 0.5, 0.5 * (4/3), (1/3) * 0.5);
@zunit = (0, 0, 0.5, 0.5);

$Lx = $La;
$Ly = $La * sqrt(3);
$Lz = $Lc;

$idx = 0;
for ($x = 0; $x < $nx; $x++)
{
  for ($y = 0; $y < $ny; $y++)
  {
    for ($z = 0; $z < $nz; $z++)
    {
      for ($i = 0; $i <= $#xunit; $i++)
      {
        $xnuc[$idx] = ($x + .25) * $Lx + $xunit[$i] * $Lx;
        $ynuc[$idx] = ($y + .25) * $Ly + $yunit[$i] * $Ly;
        $znuc[$idx] = ($z + .25) * $Lz + $zunit[$i] * $Lz;
        $idx++;
      } 
    }
  }
} 

$numnuc = $idx;

# Print length of supercell

printf("Created by Be-solid.pl\n\n");
printf("%d atoms\n",$numnuc+$numnuc*2+$numnuc*2);
printf("2 atom types\n\n");
printf("%f %f xlo xhi\n", 0, $Lx * $nx);
printf("%f %f ylo yhi\n", 0, $Ly * $ny);
printf("%f %f zlo zhi\n\n", 0, $Lz * $nz);
printf("Masses\n\n");
printf("1 9.012182\n");
printf("2 1.000000\n\n");
printf("Atoms\n\n");

$j=0;
# Print out the nuclei
for ($i = 0; $i < $numnuc; $i++)
{
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 1, 4.0, 0, 0.0, $xnuc[$i], $ynuc[$i], $znuc[$i]);
}

# Print out the core electrons
$r_core = 0.5;
for ($i = 0; $i < $numnuc; $i++)
{
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2,  0.0, 1, $r_core,$xnuc[$i], $ynuc[$i], $znuc[$i]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, -1, $r_core,$xnuc[$i], $ynuc[$i], $znuc[$i]);
}

# Print out the valence electrons
$r_valence = 2.5;
for ($i = 0; $i < $numnuc; $i += 4)
{
  $x = &BoundX($xnuc[$i] + $Lx / 2.0);
  $y = &BoundY($ynuc[$i] + $Lx / (2.0 * sqrt(3)));
  $z = $znuc[$i];
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, 1, $r_valence,$x, $y, $z);
  $x = &BoundX($xnuc[$i] + $Lx / 2.0);
  $y = &BoundY($ynuc[$i] + $Lx / (2.0 * sqrt(3)));
  $z = $znuc[$i];
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, -1, $r_valence,$x, $y, $z);

  $x = &BoundX($xnuc[$i+1] + $Lx / 2.0);
  $y = &BoundY($ynuc[$i+1] + $Lx / (2.0 * sqrt(3)));
  $z = $znuc[$i+1];
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, 1, $r_valence, $x, $y, $z);
  $x = &BoundX($xnuc[$i+1] + $Lx / 2.0);
  $y = &BoundY($ynuc[$i+1] + $Lx / (2.0 * sqrt(3)));
  $z = $znuc[$i+1];
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, -1, $r_valence,$x, $y, $z);

  $x = &BoundX($xnuc[$i+2] - $Lx / 2.0);
  $y = &BoundY($ynuc[$i+2] - $Lx / (2.0 * sqrt(3)));
  $z = $znuc[$i+2];
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, 1, $r_valence,$x, $y, $z);
  $x = &BoundX($xnuc[$i+2] - $Lx / 2.0);
  $y = &BoundY($ynuc[$i+2] - $Lx / (2.0 * sqrt(3)));
  $z = $znuc[$i+2];
  printf("%i %i %f %f %f %f %f %f\n", $j+=1, 2, 0.0, -1, $r_valence,$x, $y, $z);

  $x = &BoundX($xnuc[$i+3] - $Lx / 2.0);
  $y = &BoundY($ynuc[$i+3] - $Lx / (2.0 * sqrt(3)));
  $z = $znuc[$i+3];
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, 1, $r_valence, $x, $y, $z);
  $x = &BoundX($xnuc[$i+3] - $Lx / 2.0);
  $y = &BoundY($ynuc[$i+3] - $Lx / (2.0 * sqrt(3)));
  $z = $znuc[$i+3];
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, -1, $r_valence,$x, $y, $z);
}

sub BoundX {
  $x = $_[0];
  if ($x > $Lx * $nx) {$x -= $Lx * $nx;}
  if ($x < 0) {$x += $Lx * $nx;}
  return $x;
}

sub BoundY {
  $y = $_[0];
  if ($y > $Ly * $ny) {$y -= $Ly * $ny;}
  if ($y < 0) {$y += $Ly * $ny;}
  return $y;
}
  
sub BoundZ {
  $z = $_[0];
  if ($z > $Lz * $nz) {$z -= $Lz * $nz;}
  if ($z < 0) {$z += $Lz * $nz;}
  return $z;
}
