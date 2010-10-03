#!/usr/bin/perl

# Usage: diamond <nx> <ny> <nz> > data.diamond 
#
# Generates diamond A4 structure.

$nx = shift(@ARGV);
$ny = shift(@ARGV);
$nz = shift(@ARGV);
#$scale = shift(@ARGV);

$scale = 1.045;
$L = 7.33461;       # eFF optimized (6.740 expt)
$re_core  = 0.328;        # core electron
$re_sigma = 1.559;        # sigma electrons

#$r_threshold = $scale * 1.6 / 0.529;          # threshold for sigma bond length
$r_threshold = $scale * 1.6 / 0.5 ;          # threshold for sigma bond length

@xunit = (0, 0, 0.5, 0.5, 0.25, 0.75, 0.25, 0.75);
@yunit = (0, 0.5, 0, 0.5, 0.25, 0.75, 0.75, 0.25);
@zunit = (0, 0.5, 0.5, 0, 0.25, 0.25, 0.75, 0.75);

$idx = 0;
for ($x = 0; $x < $nx; $x++)
{
  for ($y = 0; $y < $ny; $y++)
  {
    for ($z = 0; $z < $nz; $z++)
    {
      for ($i = 0; $i <= $#xunit; $i++)
      {
        $xnuc[$idx] = $x * $L + $xunit[$i] * $L;
        $ynuc[$idx] = $y * $L + $yunit[$i] * $L;
        $znuc[$idx] = $z * $L + $zunit[$i] * $L;
        $idx++;
      } 
    }
  }
} 

$numnuc = $idx;
$xx = $L * $nx;
$yy = $L * $ny;
$zz = $L * $nz;

# Print length of supercell

printf("Created with Diamond.pl\n\n");
printf("%d atoms\n",$numnuc*7);
printf("2 atom types\n\n");
printf("%f %f xlo xhi\n", 0, $L * $nx);
printf("%f %f ylo yhi\n", 0, $L * $ny);
printf("%f %f zlo zhi\n\n", 0, $L * $nz);
printf("Masses\n\n");
printf("1 12.01070\n");
printf("2 1.000000\n\n");
printf("Atoms\n\n");

# Print out the nuclei and the core electrons
for ($i = 0; $i < $numnuc; $i++)
{
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 1, 6.0, 0.0, 0.0, $xnuc[$i], $ynuc[$i], $znuc[$i]);
}

for ($i = 0; $i < $numnuc; $i++)
{
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0, 1, $re_core,$xnuc[$i], $ynuc[$i], $znuc[$i]);
  printf("%i %i %f %i %f %f %f %f\n", $j+=1, 2, 0.0,-1, $re_core,$xnuc[$i], $ynuc[$i], $znuc[$i]);
}

# Print out sigma electrons

$LLx = $L * $nx;
$LLy = $L * $ny;
$LLz = $L * $nz;

$k=$j;
for ($i = 0; $i < $numnuc; $i++)
{
  for ($j = 0; $j < $i; $j++)
  {
    $x = $xnuc[$j] - $xnuc[$i];
    $y = $ynuc[$j] - $ynuc[$i];
    $z = $znuc[$j] - $znuc[$i];

    # Minimum image convention

    if ($x > $LLx/2) {$x -= $LLx;}
    if ($y > $LLy/2) {$y -= $LLy;}
    if ($z > $LLz/2) {$z -= $LLz;}
    
    if ($x < -$LLx/2) {$x += $LLx;}
    if ($y < -$LLy/2) {$y += $LLy;}
    if ($z < -$LLz/2) {$z += $LLz;}

    $r = sqrt($x * $x + $y * $y + $z * $z);
    if ($r   < $r_threshold)
    {
      $bond_x = $xnuc[$i] + $x / 2;
      $bond_y = $ynuc[$i] + $y / 2;
      $bond_z = $znuc[$i] + $z / 2;

      # Minimum image convention

      if ($bond_x > $LLx) {$bond_x -= $LLx};
      if ($bond_y > $LLy) {$bond_y -= $LLy};
      if ($bond_z > $LLz) {$bond_z -= $LLz};

      if ($bond_x < 0)   {$bond_x += $LLx};
      if ($bond_y < 0)   {$bond_y += $LLy};
      if ($bond_z < 0)   {$bond_z += $LLz};

      printf("%i %i %f %i %f %f %f %f\n", $k+=1, 2, 0.0, 1, $re_sigma,$bond_x, $bond_y, $bond_z);
      printf("%i %i %f %i %f %f %f %f\n", $k+=1, 2, 0.0,-1, $re_sigma,$bond_x, $bond_y, $bond_z);
    }
  }
}

