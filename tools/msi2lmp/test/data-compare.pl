#!/usr/bin/perl -w
# Tool to validate and compare two LAMMPS data files
# with "inexact" floating point comparisons
# July 2013 by Axel Kohlmeyer <akohlmey@gmail.com>

use strict;
use warnings;

my $version = 'v0.2';

# delta for floating point comparisons.
my $small = 1.0e-4;
# two hashes for storing system information
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
          die "Multiple bondangle coefficient assignments to the same angle type"
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
      } elsif ($1 eq "AngleAngleTorsion Coeffs") {
        $data->{angleangletorsioncoeff} = [];
        $i = 0;
        section_check($fh,$data,"AngleAngleTorsion Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "AngleAngleTorsion type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{ndihedraltypes}));
            ++$i;
            $data->{angleangletorsioncoeff}[$j] = $2;
            next;
          }

          die "Too many entries in AngleAngleTorsion Coeffs section"
            if ($i > $data->{ndihedraltypes});
          die "Too few entries in AngleAngleTorsion Coeffs section"
            if ($i < $data->{ndihedraltypes});
          die "Multiple dihedral coefficient assignments to the same dihedral type"
            if (scalar @{$data->{angleangletorsioncoeff}} != $data->{ndihedraltypes});

          last;
        }
      } elsif ($1 eq "EndBondTorsion Coeffs") {
        $data->{endbondtorsioncoeff} = [];
        $i = 0;
        section_check($fh,$data,"EndBondTorsion Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "EndBondTorsion type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{ndihedraltypes}));
            ++$i;
            $data->{endbondtorsioncoeff}[$j] = $2;
            next;
          }

          die "Too many entries in EndBondTorsion Coeffs section"
            if ($i > $data->{ndihedraltypes});
          die "Too few entries in EndBondTorsion Coeffs section"
            if ($i < $data->{ndihedraltypes});
          die "Multiple dihedral coefficient assignments to the same dihedral type"
            if (scalar @{$data->{endbondtorsioncoeff}} != $data->{ndihedraltypes});

          last;
        }
      } elsif ($1 eq "MiddleBondTorsion Coeffs") {
        $data->{middlebondtorsioncoeff} = [];
        $i = 0;
        section_check($fh,$data,"MiddleBondTorsion Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "MiddleBondTorsion type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{ndihedraltypes}));
            ++$i;
            $data->{middlebondtorsioncoeff}[$j] = $2;
            next;
          }

          die "Too many entries in MiddleBondTorsion Coeffs section"
            if ($i > $data->{ndihedraltypes});
          die "Too few entries in MiddleBondTorsion Coeffs section"
            if ($i < $data->{ndihedraltypes});
          die "Multiple dihedral coefficient assignments to the same dihedral type"
            if (scalar @{$data->{middlebondtorsioncoeff}} != $data->{ndihedraltypes});

          last;
        }
      } elsif ($1 eq "BondBond13 Coeffs") {
        $data->{bondbond13coeff} = [];
        $i = 0;
        section_check($fh,$data,"BondBond13 Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "BondBond13 type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{ndihedraltypes}));
            ++$i;
            $data->{bondbond13coeff}[$j] = $2;
            next;
          }

          die "Too many entries in BondBond13 Coeffs section"
            if ($i > $data->{ndihedraltypes});
          die "Too few entries in BondBond13 Coeffs section"
            if ($i < $data->{ndihedraltypes});
          die "Multiple dihedral coefficient assignments to the same dihedral type"
            if (scalar @{$data->{bondbond13coeff}} != $data->{ndihedraltypes});

          last;
        }
      } elsif ($1 eq "AngleTorsion Coeffs") {
        $data->{angletorsioncoeff} = [];
        $i = 0;
        section_check($fh,$data,"AngleTorsion Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "AngleTorsion type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{ndihedraltypes}));
            ++$i;
            $data->{angletorsioncoeff}[$j] = $2;
            next;
          }

          die "Too many entries in AngleTorsion Coeffs section"
            if ($i > $data->{ndihedraltypes});
          die "Too few entries in AngleTorsion Coeffs section"
            if ($i < $data->{ndihedraltypes});
          die "Multiple dihedral coefficient assignments to the same dihedral type"
            if (scalar @{$data->{angletorsioncoeff}} != $data->{ndihedraltypes});

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
      } elsif ($1 eq "AngleAngle Coeffs") {
        $data->{angleanglecoeff} = [];
        $i = 0;
        section_check($fh,$data,"AngleAngle Coeffs");

        while (get_next($fh)) {

          if (/^\s*([0-9]+)\s+([-+.eE0-9 ]+)\s*$/) {
            $j = $1 - 1;
            die "AngleAngle type $1 is out of range" 
              if (($1 < 1) || ($1 > $data->{nimpropertypes}));
            ++$i;
            $data->{angleanglecoeff}[$j] = $2;
            next;
          }

          die "Too many entries in AngleAngle Coeffs section"
            if ($i > $data->{nimpropertypes});
          die "Too few entries in AngleAngle Coeffs section"
            if ($i < $data->{nimpropertypes});
          die "Multiple angleangle coefficient assignments to the same angle type"
            if (scalar @{$data->{angleanglecoeff}} != $data->{nimpropertypes});

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
          if (/^\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([-+.eE0-9]+)\s+([-+.eE0-9]+)\s+([-+.eE0-9]+)\s+([-+.eE0-9]+)(|\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+))\s*$/) {

            $k = $1 - 1;
            die "Atom id $1 is out of range"
              if (($1 < 1) || ($1 > $data->{natoms}));

            $j = $3 - 1;
            die "Atom type $2 is out of range"
              if (($3 < 1) || ($3 > $data->{natomtypes}));

            ++$i;
            $data->{tag}[$k] = $1;
            $data->{molid}[$k] = $2;
            $data->{type}[$k] = $3;
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
#          } else {
#            print "Atoms: $_\n";
          }

          die "Too many entries in Atoms section: $i vs. $data->{natoms}"
            if ($i > $data->{natoms});
          die "Too few entries in Atoms section: $i vs. $data->{natoms}"
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
          if (/^\s*([0-9]+)\s+(-?[0-9]+)\s+([0-9]+)\s+([0-9]+)\s*$/) {

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
          if (/^\s*([0-9]+)\s+(-?[0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s*$/) {

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
      } elsif ($1 eq "Dihedrals") {
        $data->{dihedralt} = [];
        $data->{dihedral1} = [];
        $data->{dihedral2} = [];
        $data->{dihedral3} = [];
        $data->{dihedral4} = [];
        $i = 0;
        section_check($fh,$data,"Dihedrals");

        while (get_next($fh)) {
          if (/^\s*([0-9]+)\s+(-?[0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s*$/) {

            $k = $1 - 1;
            die "Dihedral id $1 is out of range"
              if (($1 < 1) || ($1 > $data->{ndihedrals}));

            die "Dihedral type $2 is out of range"
              if (($2 == 0) || ($2 > $data->{ndihedraltypes}));

            die "Dihedral atom 1 ID $3 is out of range"
              if (($3 < 1) || ($3 > $data->{natoms}));
            die "Dihedral atom 2 ID $4 is out of range"
              if (($4 < 1) || ($4 > $data->{natoms}));
            die "Dihedral atom 3 ID $5 is out of range"
              if (($5 < 1) || ($5 > $data->{natoms}));
            die "Dihedral atom 4 ID $6 is out of range"
              if (($6 < 1) || ($6 > $data->{natoms}));

            ++$i;
            $data->{dihedralt}[$k] = $2;
            $data->{dihedral1}[$k] = $3;
            $data->{dihedral2}[$k] = $4;
            $data->{dihedral3}[$k] = $5;
            $data->{dihedral4}[$k] = $6;
            next;
          }

          die "Too many entries in Dihedrals section"
            if ($i > $data->{ndihedrals});
          die "Too few entries in Dihedrals section"
            if ($i < $data->{ndihedrals});
          die "Multiple dihedrals assigned to the same dihedral ID"
            if (scalar @{$data->{dihedralt}} != $data->{ndihedrals});

          last;
        }
      } elsif ($1 eq "Impropers") {
        $data->{impropert} = [];
        $data->{improper1} = [];
        $data->{improper2} = [];
        $data->{improper3} = [];
        $data->{improper4} = [];
        $i = 0;
        section_check($fh,$data,"Impropers");

        while (get_next($fh)) {
          if (/^\s*([0-9]+)\s+(-?[0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s*$/) {

            $k = $1 - 1;
            die "Improper id $1 is out of range"
              if (($1 < 1) || ($1 > $data->{nimpropers}));

            die "Improper type $2 is out of range"
              if (($2 == 0) || ($2 > $data->{nimpropertypes}));

            die "Improper atom 1 ID $3 is out of range"
              if (($3 < 1) || ($3 > $data->{natoms}));
            die "Improper atom 2 ID $4 is out of range"
              if (($4 < 1) || ($4 > $data->{natoms}));
            die "Improper atom 3 ID $5 is out of range"
              if (($5 < 1) || ($5 > $data->{natoms}));
            die "Improper atom 4 ID $6 is out of range"
              if (($6 < 1) || ($6 > $data->{natoms}));

            ++$i;
            $data->{impropert}[$k] = $2;
            $data->{improper1}[$k] = $3;
            $data->{improper2}[$k] = $4;
            $data->{improper3}[$k] = $5;
            $data->{improper4}[$k] = $6;
            next;
          }

          die "Too many entries in Impropers section"
            if ($i > $data->{nimpropers});
          die "Too few entries in Impropers section"
            if ($i < $data->{nimpropers});
          die "Multiple impropers assigned to the same improper ID"
            if (scalar @{$data->{impropert}} != $data->{nimpropers});

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
  my ($n1,$n2,$rel) = @_;

  my $diff = abs($n1-$n2);
  my $avg = (abs($n1)+abs($n2))*0.5;
  return 0 if ($avg == 0.0);
  if ($rel) {
#    print "relative difference: ",$diff/$avg," vs. $small\n";
    return 0 if ($diff/$avg < $small);
  } else {
#    print "absolute difference: ",$diff," vs. $small\n";
    return 0 if ($diff < $small);
  }
  return 1;
}

sub coeffcompare {
  my ($d1,$d2,$coeff,$type) = @_;
  my (@c1,@c2,$a,$b);
  my ($field,$count,$i,$j,$t) = ($coeff . 'coeff', 'n' . $type . 'types', 0,0,0);

  if (exists $d1->{$field} && exists $d2->{$field}) {
    for ($i=0; $i < $d1->{$count}; ++$i) {
      $t = $i+1;
      @c1 = split /\s+/, ${$$d1{$field}}[$i];
      @c2 = split /\s+/, ${$$d2{$field}}[$i];
      die "Inconsistent number of $coeff coefficients for $type type $t: $#c1 vs $#c2\n"
        if ($#c1 != $#c2);

      for ($j = 0; $j <= $#c1; ++$j) {
        $a = $c1[$j]; $b = $c2[$j];
        die "Inconsistent $coeff coefficient ", $j+1,
          " for $type type $t: $a vs. $b" if (floatdiff($a,$b,1));
      }
    }
  } else {
    die "Field $field only exists in data file 1" if (exists $d1->{$field});
    die "Field $field only exists in data file 2" if (exists $d2->{$field});
  }
}

sub topocompare {
  my ($d1,$d2,$type,$count) = @_;
  my ($num,$a,$b,$field,$i,$j,$t);

  $field = 'n' . $type . 's';
  $num = $d1->{$field};

  for ($i=0; $i < $num; ++$i) {
    $t = $i+1;
    $field = $type . 't';
    $a = $d1->{$field}[$i]; $b = $d2->{$field}[$i];
    die "Inconsistent $type types for $type $t: $a vs. $b" if ($a != $b);
    for ($j=1; $j <= $count; ++$j) {
      $field = $type . $j;
      $a = $d1->{$field}[$i]; $b = $d2->{$field}[$i];
      die "Inconsistent $type atom $j for $type $t: $a vs. $b" if ($a != $b);
    }
  }
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
    die "Box data for $i does not match: $a $b" if (floatdiff($a,$b,0));
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
      die "Inconsistent data for $t, atom $j: $a vs. $b" if (floatdiff($a,$b,0));
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

  topocompare($d1,$d2,'bond',2);
  topocompare($d1,$d2,'angle',3);
  topocompare($d1,$d2,'dihedral',4);
  topocompare($d1,$d2,'improper',4);

  coeffcompare($d1,$d2,'pair','atom');
  coeffcompare($d1,$d2,'bond','bond');
  coeffcompare($d1,$d2,'angle','angle');
  coeffcompare($d1,$d2,'bondbond','angle');
  coeffcompare($d1,$d2,'bondangle','angle');
  coeffcompare($d1,$d2,'dihedral','dihedral');
  coeffcompare($d1,$d2,'angleangletorsion','dihedral');
  coeffcompare($d1,$d2,'bondbond13','dihedral');
  coeffcompare($d1,$d2,'endbondtorsion','dihedral');
  coeffcompare($d1,$d2,'middlebondtorsion','dihedral');
  coeffcompare($d1,$d2,'improper','improper');
  coeffcompare($d1,$d2,'angleangle','improper');

}

########################################################################
# main program

my $fp;

if ($#ARGV < 1) {
  die "usage $0 <file 1> <file 2>";
}

print "\nLAMMPS data file validation tool. $version\n\n";

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
