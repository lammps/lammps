#!/usr/bin/perl -w

use strict;
use warnings;

# delta for floating point comparisons.
my $small = 1.0e-10;
my $invalid = -9.8765321e100;
# hashes for system information
my %data1;
my %data2;

# simple checks after reading a section header
sub section_check {
  my ($fh,$data,$section) = @_;
  $_ = <$fh>;
  # skip empty and whitespace-only lines
  die "Line following '$section' is not empty" unless (/^\s*$/);
  die "Incomplete or incorrect header" unless ($data->{natoms} > 0);
  die "Incomplete or incorrect header" unless ($data->{natomtypes} > 0);
}

sub get_next {
  my ($fh) = @_;

  while (<$fh>) {
    chomp;
    # trim off comments
    $_ =~ s/#.*$//;
    # skip empty and whitespace-only lines
    next if (/^\s*$/);
    last;
  }
  return $_;
}

# fill hash with default data.
sub data_defaults {
  my ($data) = @_;

  $data->{natoms} = 0;
  $data->{nbonds} = 0;
  $data->{nangles} = 0;
  $data->{ndihedrals} = 0;
  $data->{nimpropers} = 0;

  $data->{natomtypes} = 0;
  $data->{nbondtypes} = 0;
  $data->{nangletypes} = 0;
  $data->{ndihedraltypes} = 0;
  $data->{nimpropertypes} = 0;

  $data->{xlo} =  0.5;
  $data->{xhi} = -0.5;
  $data->{ylo} =  0.5;
  $data->{yhi} = -0.5;
  $data->{zlo} =  0.5;
  $data->{zhi} = -0.5;
  $data->{xy} =   0.0;
  $data->{xz} =   0.0;
  $data->{yz} =   0.0;
  $data->{triclinic} = 0;
}

# read/parse lammps data file
sub read_data {
  my ($fh,$data) = @_;
  my $section;

  # read header. first line is already chopped off
  while (get_next($fh)) {

    if (/^\s*([0-9]+)\s+atoms\s*$/)     { $data->{natoms} = $1;     next; }
    if (/^\s*([0-9]+)\s+bonds\s*$/)     { $data->{nbonds} = $1;     next; }
    if (/^\s*([0-9]+)\s+angles\s*$/)    { $data->{nangles} = $1;    next; }
    if (/^\s*([0-9]+)\s+dihedrals\s*$/) { $data->{ndihedrals} = $1; next; }
    if (/^\s*([0-9]+)\s+impropers\s*$/) { $data->{nimpropers} = $1; next; }

    if (/^\s*([0-9]+)\s+atom types\s*$/)   { $data->{natomtypes} = $1;  next; }
    if (/^\s*([0-9]+)\s+bond types\s*$/)   { $data->{nbondtypes} = $1;  next; }
    if (/^\s*([0-9]+)\s+angle types\s*$/)  { $data->{nangletypes} = $1; next; }
    if (/^\s*([0-9]+)\s+dihedral types\s*$/)
      { $data->{ndihedraltypes} = $1; next; }
    if (/^\s*([0-9]+)\s+improper types\s*$/)
      { $data->{nimpropertypes} =$1 ; next; }

    if (/^\s*([-+.eE0-9]+)\s+([-+.eE0-9]+)\s+xlo xhi\s*$/) {
      $data->{xlo}=$1;
      $data->{xhi}=$2;
      next;
    }
    if (/^\s*([-+.eE0-9]+)\s+([-+.eE0-9]+)\s+ylo yhi\s*$/) {
      $data->{ylo}=$1;
      $data->{yhi}=$2;
      next;
    }
    if (/^\s*([-+.eE0-9]+)\s+([-+.eE0-9]+)\s+zlo zhi\s*$/) {
      $data->{zlo}=$1;
      $data->{zhi}=$2;
      next;
    }
    if (/^\s*([-+.eE0-9]+)\s+([-+.eE0-9]+)\s+([-+.eE0-9]+)\s+xy xz yz\s*$/) {
      $data->{xy}=$1;
      $data->{xz}=$2;
      $data->{yz}=$3;
      $data->{triclinic} = 1;
      next;
    }

    # if we reach this point, the header has ended;
    last;
  }

  $a = $data->{natoms};
  $b = $data->{natomtypes};
  die "Invalid number of atoms: $a" unless ($a > 0);
  die "Invalid number of atom types: $b" unless ($b > 0);

  my ($i,$j,$k);
  while (1) {
    if (/^\s*(\S+|\S+ Coeffs)\s*$/) {
      if ($1 eq "Masses") {
        $data->{mass} = [];
        $i = 0;
        section_check($fh,$data,"Masses");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9]+)\s*$/) {
            $j = $1 - 1;
            die "Atom type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{natomtypes}));
            ++$i;
            $data->{mass}[$j] = $2;
            next;
          }

          die "Too many entries in Masses section"
            if ($i > $data->{natomtypes});
          die "Too few entries in Masses section"
            if ($i < $data->{natomtypes});
          die "Multiple mass assignments to the same atom type"
            if (scalar @{$data->{mass}} != $data->{natomtypes});

          last;
        }
      } elsif ($1 eq "Pair Coeffs") {
        $data->{paircoeff} = [];
        $i = 0;
        section_check($fh,$data,"Pair Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "Atom type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{natomtypes}));
            ++$i;
            $data->{paircoeff}[$j] = $2;
            next;
          }

          die "Too many entries in Pair Coeffs section"
            if ($i > $data->{natomtypes});
          die "Too few entries in Pair Coeffs section"
            if ($i < $data->{natomtypes});
          die "Multiple pair coefficient assignments to the same atom type"
            if (scalar @{$data->{paircoeff}} != $data->{natomtypes});

          last;
        }
      } elsif ($1 eq "Bond Coeffs") {
        $data->{bondcoeff} = [];
        $i = 0;
        section_check($fh,$data,"Bond Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "Bond type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{nbondtypes}));
            ++$i;
            $data->{bondcoeff}[$j] = $2;
            next;
          }

          die "Too many entries in Bond Coeffs section"
            if ($i > $data->{nbondtypes});
          die "Too few entries in Bond Coeffs section"
            if ($i < $data->{nbondtypes});
          die "Multiple bond coefficient assignments to the same bond type"
            if (scalar @{$data->{bondcoeff}} != $data->{nbondtypes});

          last;
        }
      } elsif ($1 eq "Angle Coeffs") {
        $data->{anglecoeff} = [];
        $i = 0;
        section_check($fh,$data,"Angle Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "Angle type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{nangletypes}));
            ++$i;
            $data->{anglecoeff}[$j] = $2;
            next;
          }

          die "Too many entries in Angle Coeffs section"
            if ($i > $data->{nangletypes});
          die "Too few entries in Angle Coeffs section"
            if ($i < $data->{nangletypes});
          die "Multiple angle coefficient assignments to the same angle type"
            if (scalar @{$data->{anglecoeff}} != $data->{nangletypes});

          last;
        }
      } elsif ($1 eq "BondBond Coeffs") {
        $data->{bondbondcoeff} = [];
        $i = 0;
        section_check($fh,$data,"BondBond Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "Angle type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{nangletypes}));
            ++$i;
            $data->{bondbondcoeff}[$j] = $2;
            next;
          }

          die "Too many entries in BondBond Coeffs section"
            if ($i > $data->{nangletypes});
          die "Too few entries in BondBond Coeffs section"
            if ($i < $data->{nangletypes});
          die "Multiple angle coefficient assignments to the same angle type"
            if (scalar @{$data->{bondbondcoeff}} != $data->{nangletypes});

          last;
        }
      } elsif ($1 eq "BondAngle Coeffs") {
        $data->{bondanglecoeff} = [];
        $i = 0;
        section_check($fh,$data,"BondAngle Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "Angle type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{nangletypes}));
            ++$i;
            $data->{bondanglecoeff}[$j] = $2;
            next;
          }

          die "Too many entries in BondAngle Coeffs section"
            if ($i > $data->{nangletypes});
          die "Too few entries in BondAngle Coeffs section"
            if ($i < $data->{nangletypes});
          die "Multiple angle coefficient assignments to the same angle type"
            if (scalar @{$data->{bondanglecoeff}} != $data->{nangletypes});

          last;
        }
      } elsif ($1 eq "Dihedral Coeffs") {
        $data->{dihedralcoeff} = [];
        $i = 0;
        section_check($fh,$data,"Dihedral Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "Dihedral type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{ndihedraltypes}));
            ++$i;
            $data->{dihedralcoeff}[$j] = $2;
            next;
          }

          die "Too many entries in Dihedral Coeffs section"
            if ($i > $data->{ndihedraltypes});
          die "Too few entries in Dihedral Coeffs section"
            if ($i < $data->{ndihedraltypes});
          die "Multiple dihedral coefficient assignments to the same dihedral type"
            if (scalar @{$data->{dihedralcoeff}} != $data->{ndihedraltypes});

          last;
        }
      } elsif ($1 eq "Improper Coeffs") {
        $data->{impropercoeff} = [];
        $i = 0;
        section_check($fh,$data,"Improper Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "Improper type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{nimpropertypes}));
            ++$i;
            $data->{impropercoeff}[$j] = $2;
            next;
          }

          die "Too many entries in Improper Coeffs section"
            if ($i > $data->{nimpropertypes});
          die "Too few entries in Improper Coeffs section"
            if ($i < $data->{nimpropertypes});
          die "Multiple improper coefficient assignments to the same improper type"
            if (scalar @{$data->{impropercoeff}} != $data->{nimpropertypes});

          last;
        }
      } elsif ($1 eq "Atoms") {
        $data->{tag} = [];
        $data->{type} = [];
        $data->{molid} = [];
        $data->{charge} = [];
        $data->{posx} = [];
        $data->{posy} = [];
        $data->{posz} = [];
        $data->{imgx} = [];
        $data->{imgy} = [];
        $data->{imgz} = [];
        $i = 0;
        section_check($fh,$data,"Atoms");

        while (get_next($fh)) {
          if (/^\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([-+.eE0-9]+)\s+([-+.eE0-9]+)\s+([-+.eE0-9]+)\s+([-+.eE0-9]+)(|\s+([0-9]+)\s+([0-9]+)\s+([0-9]+))\s*$/) {

            $k = $1 - 1;
            die "Atom id $1 is out of range"
              if (($1 < 1) || ($1 > $data->{natoms}));

            $j = $2 - 1;
            die "Atom type $2 is out of range"
              if (($2 < 1) || ($2 > $data->{natomtypes}));

            ++$i;
            $data->{tag}[$k] = $1;
            $data->{type}[$k] = $2;
            $data->{molid}[$k] = $3;
            $data->{charge}[$k] = $4;
            $data->{posx}[$k] = $5;
            $data->{posy}[$k] = $6;
            $data->{posz}[$k] = $7;
            $data->{imgx}[$k] = 0;
            $data->{imgy}[$k] = 0;
            $data->{imgz}[$k] = 0;
            if (! $8 eq "") {
              $data->{imgx}[$k] = $9;
              $data->{imgy}[$k] = $10;
              $data->{imgz}[$k] = $11;
            }
            next;
          }

          die "Too many entries in Atoms section"
            if ($i > $data->{natoms});
          die "Too few entries in Atoms section"
            if ($i < $data->{natoms});
          die "Multiple atoms assigned to the same atom ID"
            if (scalar @{$data->{tag}} != $data->{natoms});

          last;
        }
      } elsif ($1 eq "Velocities") {
        $data->{velx} = [];
        $data->{vely} = [];
        $data->{velz} = [];
        $i = 0;
        section_check($fh,$data,"Velocities");

        while (get_next($fh)) {
          if (/^\s*([0-9]+)\s+([-+.eE0-9]+)\s+([-+.eE0-9]+)\s+([-+.eE0-9]+)\s*$/) {

            $k = $1 - 1;
            die "Atom id $1 is out of range"
              if (($1 < 1) || ($1 > $data->{natoms}));

            ++$i;
            $data->{velx}[$k] = $2;
            $data->{vely}[$k] = $3;
            $data->{velz}[$k] = $4;
            next;
          }

          die "Too many entries in Velocities section"
            if ($i > $data->{natoms});
          die "Too few entries in Velocities section"
            if ($i < $data->{natoms});
          die "Multiple velocities assigned to the same atom ID"
            if (scalar @{$data->{velx}} != $data->{natoms});

          last;
        }
      } elsif ($1 eq "Bonds") {
        $data->{bondt} = [];
        $data->{bond1} = [];
        $data->{bond2} = [];
        $i = 0;
        section_check($fh,$data,"Bonds");

        while (get_next($fh)) {
          if (/^\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s*$/) {

            $k = $1 - 1;
            die "Bond id $1 is out of range"
              if (($1 < 1) || ($1 > $data->{nbonds}));

            die "Bond type $2 is out of range"
              if (($2 == 0) || ($2 > $data->{nbondtypes}));

            die "Bond atom 1 ID $3 is out of range"
              if (($3 < 1) || ($3 > $data->{natoms}));
            die "Bond atom 2 ID $4 is out of range"
              if (($4 < 1) || ($4 > $data->{natoms}));

            ++$i;
            $data->{bondt}[$k] = $2;
            $data->{bond1}[$k] = $3;
            $data->{bond2}[$k] = $4;
            next;
          }

          die "Too many entries in Bonds section"
            if ($i > $data->{nbonds});
          die "Too few entries in Bonds section"
            if ($i < $data->{nbonds});
          die "Multiple bonds assigned to the same bond ID"
            if (scalar @{$data->{bondt}} != $data->{nbonds});

          last;
        }
      } elsif ($1 eq "Angles") {
        $data->{anglet} = [];
        $data->{angle1} = [];
        $data->{angle2} = [];
        $data->{angle3} = [];
        $i = 0;
        section_check($fh,$data,"Angles");

        while (get_next($fh)) {
          if (/^\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s*$/) {

            $k = $1 - 1;
            die "Angle id $1 is out of range"
              if (($1 < 1) || ($1 > $data->{nangles}));

            die "Angle type $2 is out of range"
              if (($2 == 0) || ($2 > $data->{nangletypes}));

            die "Angle atom 1 ID $3 is out of range"
              if (($3 < 1) || ($3 > $data->{natoms}));
            die "Angle atom 2 ID $4 is out of range"
              if (($4 < 1) || ($4 > $data->{natoms}));
            die "Angle atom 3 ID $5 is out of range"
              if (($5 < 1) || ($5 > $data->{natoms}));

            ++$i;
            $data->{anglet}[$k] = $2;
            $data->{angle1}[$k] = $3;
            $data->{angle2}[$k] = $4;
            $data->{angle3}[$k] = $5;
            next;
          }

          die "Too many entries in Angles section"
            if ($i > $data->{nangles});
          die "Too few entries in Angles section"
            if ($i < $data->{nangles});
          die "Multiple angles assigned to the same angle ID"
            if (scalar @{$data->{anglet}} != $data->{nangles});

          last;
        }
      } else {
        die "Bad data: $_";
      }

      last unless ($_);
    } else {
      die "Bad data: $_";
    }
  }
}

sub floatdiff {
  my ($n1,$n2) = @_;

  my $diff = abs($n1-$n2);
  my $avg = (abs($n1)+abs($n2))/2.0;
  return 0 if ($avg == 0.0);
  return 0 if ($diff < $small);
  return 1;
}

sub syscompare {
  my ($d1,$d2) = @_;
  my ($i,$j,$t,$a,$b,@l);

  # check atoms.
  die "Number of atoms does not match"
    if ($d1->{natoms} != $d2->{natoms});
  die "Number of atom types does not match"
    if ($d1->{natomtypes} != $d2->{natomtypes});

  # check bonded interactions
  @l = ('bond','angle','dihedral','improper');
  foreach $i (@l) {
    $t = 'n' . $i . 's';
    $a = $d1->{$t};
    $b = $d2->{$t};
    die "Number of ",$i,"s does not match: $a vs $b" unless ($a == $b);

    $t = 'n' . $i . 'types';
    $a = $d1->{$t};
    $b = $d2->{$t};
    die "Number of ",$i," types does not match: $a vs $b" unless ($a == $b);
  }

  # check box information
  die "Inconsistent box shape" if ($d1->{triclinic} != $d2->{triclinic});

  @l = ('xlo','xhi','ylo','yhi','zlo','zhi');
  if ($d1->{triclinic}) { push @l, ('xy','xz','yz'); }
  for $i (@l) {
    $a = $d1->{$i};
    $b = $d2->{$i};
    die "Box data for $i does not match: $a $b" if (floatdiff($a,$b));
  }

  for ($i=0; $i < $d1->{natoms}; ++$i) {
    $j = $i+1;
    for $t ('tag','molid','type','imgx','imgy','imgz') {
      if (exists $d1->{$t}[$i]) {
        $a = $d1->{$t}[$i];
      } else {
        $a = 0;
      }
      if (exists $d2->{$t}[$i]) {
        $b = $d2->{$t}[$i];
      } else {
        $b = 0;
      }
      die "Inconsistent data for $t, atom $j: $a vs. $b" if ($a != $b);
    }

    for $t ('charge','posx','posy','posz') {
      if (exists $d1->{$t}[$i]) {
        $a = $d1->{$t}[$i];
      } else {
        $a = 0;
      }
      if (exists $d2->{$t}[$i]) {
        $b = $d2->{$t}[$i];
      } else {
        $b = 0;
      }
      die "Inconsistent data for $t, atom $j: $a vs. $b" if ($a != $b);
    }
  }

  for ($i=0; $i < $d1->{natomtypes}; ++$i) {
    $j = $i+1;
    if (exists $d1->{mass}[$i]) {
      $a = $d1->{mass}[$i];
    } else {
      die "No mass for atom type $j in data file 1";
    }
    if (exists $d2->{mass}[$i]) {
      $a = $d2->{mass}[$i];
    } else {
      die "No mass for atom type $j in data file 2";
    }
  }

  my (@c1,@c2);
  if (exists $d1->{paircoeff} && exists $d2->{paircoeff}) {
    for ($i=0; $i < $d1->{natomtypes}; ++$i) {
      $t = $i+1;
      @c1 = split /\s+/, ${$$d1{paircoeff}}[$i];
      @c2 = split /\s+/, ${$$d2{paircoeff}}[$i];
      die "Inconsistent number of pair coefficients for atom type $t: @c1 vs @c2\n"
        if ($#c1 != $#c2);

      for ($j = 0; $j <= $#c1; ++$j) {
        die "Inconsistent pair coefficient ", $j+1, " for atom type $t"
          if (floatdiff($c1[$j],$c2[$j]));
      }
    }
  }

  if (exists $d1->{bondcoeff} && exists $d2->{bondcoeff}) {
    for ($i=0; $i < $d1->{nbondtypes}; ++$i) {
      $t = $i+1;
      @c1 = split /\s+/, ${$$d1{bondcoeff}}[$i];
      @c2 = split /\s+/, ${$$d2{bondcoeff}}[$i];
      die "Inconsistent number of bond coefficients for bond type $t: @c1 vs @c2\n"
        if ($#c1 != $#c2);

      for ($j = 0; $j <= $#c1; ++$j) {
        die "Inconsistent bond coefficient ", $j+1, " for bond type $t"
          if (floatdiff($c1[$j],$c2[$j]));
      }
    }
  }

  if (exists $d1->{anglecoeff} && exists $d2->{anglecoeff}) {
    for ($i=0; $i < $d1->{nangletypes}; ++$i) {
      $t = $i+1;
      @c1 = split /\s+/, ${$$d1{anglecoeff}}[$i];
      @c2 = split /\s+/, ${$$d2{anglecoeff}}[$i];
      die "Inconsistent number of angle coefficients for angle type $t: @c1 vs @c2\n"
        if ($#c1 != $#c2);

      for ($j = 0; $j <= $#c1; ++$j) {
        die "Inconsistent angle coefficient ", $j+1, " for angle type $t"
          if (floatdiff($c1[$j],$c2[$j]));
      }
    }
  }

  if (exists $d1->{dihedralcoeff} && exists $d2->{dihedralcoeff}) {
    for ($i=0; $i < $d1->{ndihedraltypes}; ++$i) {
      $t = $i+1;
      @c1 = split /\s+/, ${$$d1{dihedralcoeff}}[$i];
      @c2 = split /\s+/, ${$$d2{dihedralcoeff}}[$i];
      die "Inconsistent number of dihedral coefficients for dihedral type $t: @c1 vs @c2\n"
        if ($#c1 != $#c2);

      for ($j = 0; $j <= $#c1; ++$j) {
        die "Inconsistent dihedral coefficient ", $j+1, " for dihedral type $t"
          if (floatdiff($c1[$j],$c2[$j]));
      }
    }
  }

  if (exists $d1->{impropercoeff} && exists $d2->{impropercoeff}) {
    for ($i=0; $i < $d1->{nimpropertypes}; ++$i) {
      $t = $i+1;
      @c1 = split /\s+/, ${$$d1{impropercoeff}}[$i];
      @c2 = split /\s+/, ${$$d2{impropercoeff}}[$i];
      die "Inconsistent number of improper coefficients for improper type $t: @c1 vs @c2\n"
        if ($#c1 != $#c2);

      for ($j = 0; $j <= $#c1; ++$j) {
        die "Inconsistent improper coefficient ", $j+1, " for improper type $t"
          if (floatdiff($c1[$j],$c2[$j]));
      }
    }
  }
}

########################################################################
# main program

my $fp;

if ($#ARGV < 1) {
  die "usage $0 <file 1> <file 2>";
}

print "\nLAMMPS data file validation tool. v0.1\n\n";

data_defaults(\%data1);
data_defaults(\%data2);

# read in first data file
open($fp, '<', $ARGV[0]) or die $!;
print "opened data file 1: $ARGV[0]\n";
$_=<$fp>;
print;
read_data($fp,\%data1);
print "done reading data file 1\n\n";
close $fp;

# read in second data file
open($fp, '<', $ARGV[1]) or die $!;
print "opened data file 2: $ARGV[1]\n";
$_=<$fp>;
print;
read_data($fp,\%data2);
print "done reading data file 2\n\n";
close $fp;

# compare data sets
syscompare(\%data1,\%data2);

print "File $ARGV[0] and $ARGV[1] match\n\n";

exit 0;
