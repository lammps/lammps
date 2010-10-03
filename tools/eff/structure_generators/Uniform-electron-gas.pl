#!/usr/bin/perl

# Usage: ./Uniform-electron-gas.pl <nx> <ny> <nz> <rs> > data.uniform-e-gas
#
# Generates a uniform electron gas using electrons on an NaCl lattice.

$nx = shift(@ARGV);
$ny = shift(@ARGV);
$nz = shift(@ARGV);
$rs = shift(@ARGV);

# This part changes for different lattices
@xunit = (0, 0.5, 0.5, 0, 0, 0.5, 0.5, 0);
@yunit = (0, 0, 0.5, 0.5, 0.5, 0.5, 0, 0);
@zunit = (0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5);
$volume = 1;
$elecs_per_cell = 8;

$L = $rs * ((4/3) * 3.14159265 * $elecs_per_cell / $volume) ** (1/3);  # length of unit cell
$r_elec = $rs * sqrt(1.5 / 1.1050);

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
        $xnuc[$idx] = $x * $Lx + $xunit[$i] * $L;
        $ynuc[$idx] = $y * $Ly + $yunit[$i] * $L;
        $znuc[$idx] = $z * $Lz + $zunit[$i] * $L;
        $idx++;
      } 
    }
  }
} 

$numnuc = $idx;

# Print length of supercell

printf("Created with Uniform-electron-gas.pl\n\n");
printf("%d atoms\n",$numnuc);
printf("1 atom types\n\n");
printf("%f %f xlo xhi\n", 0, $Lx * $nx);
printf("%f %f ylo yhi\n", 0, $Ly * $ny);
printf("%f %f zlo zhi\n\n", 0, $Lz * $nz);
printf("Masses\n\n");
printf("1 1.000000\n\n");
printf("Atoms\n\n");

$j = 0;

# Print out the electrons
for ($i = 0; $i < $numnuc; $i += 2)
{
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 0.0, 1, 1, $r_elec, $xnuc[$i], $ynuc[$i], $znuc[$i]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 0.0, 1, -1, $r_elec, $xnuc[$i+1], $ynuc[$i+1], $znuc[$i+1]);
}

printf("\n");
