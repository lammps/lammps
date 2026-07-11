#!/usr/bin/perl

# Usage: h2.pl <nx> <ny> <nz> <rs> > data.h2
#
# Generates rectangular lattice of hydrogen atoms, # in each direction = nx,
# with Wigner radius rs.
#    ...
# Shifts H atoms in alternating directions so they form H2 molecule
# starting structures.

$nx = shift(@ARGV);
$ny = shift(@ARGV);
$nz = shift(@ARGV);
$rs = shift(@ARGV);

$L  = 1.6119919540164693 * $rs;  # length of unit cell ((4/3 pi)^(1/3)
$re = 1.823572;                  # electron radius

$dshift = 0.5;

$idx = 0;
for ($x = 0; $x < $nx; $x++)
{
  for ($y = 0; $y < $ny; $y++)
  {
    $dsign = 1; 
    for ($z = 0; $z < $nz; $z++)
    {
      $xnuc[$idx] = $x * $L + 0.5 * $L; 
      $ynuc[$idx] = $y * $L + 0.5 * $L;
      $znuc[$idx] = $z * $L + 0.5 * $L + $dshift * $dsign;
      $dsign = -$dsign;
      $idx++;
    }
  }
} 

$numnuc = $idx;

# Print length of supercell

printf("Created with h2.pl\n\n");
printf("%d atoms\n",$numnuc*2);
printf("2 atom types\n\n");
printf("%f %f xlo xhi\n", 0, $L * $nx);
printf("%f %f ylo yhi\n", 0, $L * $ny);
printf("%f %f zlo zhi\n\n", 0, $L * $nz);
printf("Masses\n\n");
printf("1 1.007940\n");
printf("2 1.000000\n\n");
printf("Atoms\n\n");

$j=0;

# Print out the nuclei and the core electrons
for ($i = 0; $i < $numnuc; $i++)
{
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 1, 1.0, 0, 0.0, $xnuc[$i], $ynuc[$i], $znuc[$i]);
}

$spin = 1;
for ($i = 0; $i < $numnuc; $i++)
{
  if ($spin == 1) {$spin = -1;} else {$spin = 1;}
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2,  0.0, $spin, $re, $xnuc[$i], $ynuc[$i], $znuc[$i]);
}

