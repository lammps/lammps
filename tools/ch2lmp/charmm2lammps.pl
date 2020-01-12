#!/usr/bin/perl
#
#  program:     charmm2lammps.pl
#  author:      Pieter J. in 't Veld,
#               pjintve@sandia.gov, veld@verizon.net
#  date:        February 12-23, April 5, 2005.
#  purpose:     Translation of charmm input to lammps input
#
#  Notes:       Copyright by author for Sandia National Laboratories
#    20050212   Needed (in the same directory):
#               - $project.crd          ; Assumed to be correct and running
#               - $project.psf          ; CHARMM configs
#               - top_$forcefield.rtf   ;
#               - par_$forcefield.prm   ;
#               Ouput:
#               - $project.data         ; LAMMPS data file
#               - $project.in           ; LAMMPS input file
#               - $project_ctrl.pdb     ; PDB control file
#               - $project_ctrl.psf     ; PSF control file
#    20050218   Optimized for memory usage
#    20050221   Rotation added
#    20050222   Water added
#    20050223   Ions added
#    20050405   Water bug fixed; addition of .pdb input
#    20050407   project_ctrl.psf bug fixed; addition of -border
#    20050519   Added interpretation of charmm xplor psfs
#    20050603   Fixed centering issues
#    20050630   Fixed symbol issues arising from salt addition
#    20060818   Changed reading of pdb format to read exact columns
#    20070109   Changed AddMass() to use $max_id correctly
#    20160114   Added compatibility for parameter files that use IMPROPERS instead of IMPROPER
#               Print warning when not all parameters are detected. Set correct number of atom types.
#    20160613   Fix off-by-one issue in atom type validation check
#               Replace -charmm command line flag with -nohints flag
#               and enable type hints in data file by default.
#               Add hints also to section headers
#               Add a brief minimization to example input template.
#    20161001   Added 'CMAP crossterms' section at the end of the data file
#    20161001   Added instructions in CMAP section to fix problem if 'ter'
#                 is not designated in the .pdb file to identify last amino acid
#    20161005   Added tweak to embed command line in generated LAMMPS input
#    20181120   Fix topology parsing bug
#
#    General    Many thanks to Paul S. Crozier for checking script validity
#               against his projects.
#               Also thanks to Xiaohu Hu (hux2@ornl.gov) and Robert A. Latour
#               (latourr@clemson.edu), David Hyde-Volpe, and Tigran Abramyan,
#               Clemson University and Chris Lorenz (chris.lorenz@kcl.ac.uk),
#               King's College London for their efforts to add CMAP sections,
#               which is implemented using the option flag "-cmap".

# Initialization

  sub Test
  {
    my $name            = shift(@_);

    printf("Error: file %s not found\n", $name) if (!scalar(stat($name)));
    return !scalar(stat($name));
  }


  sub Initialize                                        # initialization
  {
    my $k               = 0;
    my @dir             = ("x", "y", "z");
    my @options         = ("-help", "-nohints", "-water", "-ions", "-center",
                           "-quiet", "-pdb_ctrl", "-l", "-lx", "-ly", "-lz",
                           "-border", "-ax", "-ay", "-az", "-cmap");
    my @remarks         = ("display this message",
                           "do not print type and style hints in data file",
                           "add TIP3P water [default: 1 g/cc]",
                           "add (counter)ions using Na+ and Cl- [default: 0 mol/l]",
                           "recenter atoms",
                           "do not print info",
                           "output project_ctrl.pdb [default: on]",
                           "set x-, y-, and z-dimensions simultaneously",
                           "x-dimension of simulation box",
                           "y-dimension of simulation box",
                           "z-dimension of simulation box",
                           "add border to all sides of simulation box [default: 0 A]",
                           "rotation around x-axis",
                           "rotation around y-axis",
                           "rotation around z-axis",
                           "generate a CMAP section in data file"
                         );
    my $notes;

    $program            = "charmm2lammps";
    $version            = "1.9.2";
    $year               = "2018";
    $add                = 1;
    $water_dens         = 0;
    $ions               = 0;
    $info               = 1;
    $center             = 0;
    $net_charge         = 0;
    $ion_molar          = 0;
    $pdb_ctrl           = 1;
    $border             = 0;
    $L                  = (0, 0, 0);
    $cmap               = 0;
    @R                  = M_Unit();

    $notes  = "  * The average of extremes is used as the origin\n";
    $notes .= "  * Residues are numbered sequentially\n";
    $notes .= "  * Water is added on an FCC lattice: allow 5 ps for";
    $notes .= " equilibration\n";
    $notes .= "  * Ions are added randomly and only when water is present\n";
    $notes .= "  * CHARMM force field v2.7 parameters used for";
    $notes .= " water and NaCl\n";
    $notes .= "  * Rotation angles are in degrees\n";
    $notes .= "  * Rotations are executed consecutively: -ax -ay != -ay -ax\n";
    $notes .= "  * CHARMM files needed in execution directory:\n";
    $notes .= "      - project.crd           coordinates\n";
    $notes .= "      - project.pdb           when project.crd is absent\n";
    $notes .= "      - project.psf           connectivity\n";
    $notes .= "      - top_forcefield.rtf    topology\n";
    $notes .= "      - par_forcefield.prm    parameters\n";
    $notes .= "  * Output files written to execution directory:\n";
    $notes .= "      - project.data          LAMMPS data file\n";
    $notes .= "      - project.in            suggested LAMMPS input script\n";
    $notes .= "      - project_ctrl.pdb      control file when requested\n";

    # record full command line for later use
    $cmdline = "$program.pl " . join(" ",@ARGV);

    foreach (@ARGV)
    {
      if (substr($_, 0, 1) eq "-")
      {
        my $k           = 0;
        my @tmp         = split("=");
        my $switch      = ($arg[1] eq "")||($arg[1] eq "on")||($arg[1]!=0);
        $tmp[0]         = lc($tmp[0]);
        foreach (@options)
        {
          last if ($tmp[0] eq substr($_, 0 , length($tmp[0])));
          ++$k;
        }
        $help           = 1 if (!$k--);                              # -help
        $add            = 0 if (!$k--);                              # -nohints
        $water_dens     = ($tmp[1] ne "" ? $tmp[1] : 1) if (!$k--);  # -water
        $ion_molar      = abs($tmp[1]) if (!$k);                     # -ions
        $ions           = 1 if (!$k--);                              # ...
        $center         = 1 if (!$k--);                              # -center
        $info           = 0 if (!$k--);                              # -quiet
        $pdb_ctrl       = $switch if (!$k--);                        # -pdb_ctrl
        my $flag        = $k--;                                      # -l
        $L[0]           = abs($tmp[1]) if (!($flag && $k--));        # -lx
        $L[1]           = abs($tmp[1]) if (!($flag && $k--));        # -ly
        $L[2]           = abs($tmp[1]) if (!($flag && $k--));        # -lz
        $border         = abs($tmp[1]) if (!$k--);                   # -border
        @R              = M_Dot(M_Rotate(0, $tmp[1]), @R) if (!$k--);# -ax
        @R              = M_Dot(M_Rotate(1, $tmp[1]), @R) if (!$k--);# -ay
        @R              = M_Dot(M_Rotate(2, $tmp[1]), @R) if (!$k--);# -az
        $cmap           = ($tmp[1] ne "" ? $tmp[1] : 22) if (!$k--); # -cmap
        print("Warning: ignoring unknown command line flag: $tmp[0]\n") unless $k;
      }
      else
      {
        $forcefield     = $_ if (!$k);
        $project        = $_ if ($k++ == 1);
      }
    }
    $water_dens         = 1 if ($ions && !$water_dens);
    if (($k<2)||$help)
    {
      printf("\n%s v%s (c)2005-%s by Pieter J. in \'t Veld and others\n\n",
        $program, $version, $year);
      printf("Usage:\n  %s.pl [-option[=#] ..] forcefield project\n\n",$program);
      printf("Options:\n");
      for (my $i=0; $i<scalar(@options); ++$i)
      {
        printf("  %-10.10s %s\n", $options[$i], $remarks[$i]);
      }
      printf("\nNotes:\n%s\n", $notes);
      exit(-1);
    }
    else { printf("\n%s v%s (c)2005-%s\n\n", $program, $version, $year) if ($info); }
    my $flag            = Test($Parameters = "par_$forcefield.prm");
    $flag               |= Test($Topology = "top_$forcefield.rtf");
    $flag               |= Test($Pdb = "$project.pdb")
                                    if (!scalar(stat($Crd = "$project.crd")));
    $flag               |= Test($Psf = "$project.psf") if ($look eq "");
    $pdb                = ($Pdb ne "") ? 1 : 0;
    printf("Conversion aborted\n\n") if ($flag);
    exit(-1) if ($flag);
    printf("Info: using $Pdb instead of $Crd\n") if (!scalar(stat($Crd)));
    for (my $i=0; $i<3; ++$i)
    {
      printf("Info: l%s not set: will use extremes\n",
        ("x", "y", "z")[$i]) if ($info&&!$L[$i]);
    }
    open(PARAMETERS, "<par_$forcefield.prm");
  }


# Vector manipulation

  sub V_String
  {
    my @v               = @_;

    return "{".$v[0].", ".$v[1].", ".$v[2]."}";
  }


  sub V_Add
  {
    my @v1              = splice(@_, 0, 3);
    my @v2              = splice(@_, 0, 3);

    return ($v1[0]+$v2[0], $v1[1]+$v2[1], $v1[2]+$v2[2]);
  }


  sub V_Subtr
  {
    my @v1              = splice(@_, 0, 3);
    my @v2              = splice(@_, 0, 3);

    return ($v1[0]-$v2[0], $v1[1]-$v2[1], $v1[2]-$v2[2]);
  }


  sub V_Dot
  {
    my @v1              = splice(@_, 0, 3);
    my @v2              = splice(@_, 0, 3);

    return $v1[0]*$v2[0]+$v1[1]*$v2[1]+$v1[2]*$v2[2];
  }


  sub V_Mult
  {
    my @v               = splice(@_, 0, 3);
    my $f               = shift(@_);

    return ($f*$v[0], $f*$v[1], $f*$v[2]);
  }


  sub M_String
  {
    my $string;

    for (my $i=0; $i<3; ++$i)
    {
      $string           .= ", " if ($i);
      $string           .= V_String(splice(@_, 0, 3));
    }
    return "{".$string."}";
  }


  sub M_Transpose
  {
    return
      (@_[0], @_[3], @_[6],
       @_[1], @_[4], @_[7],
       @_[2], @_[5], @_[8]);
  }


  sub M_Dot
  {
    my @v11             = splice(@_, 0, 3);
    my @v12             = splice(@_, 0, 3);
    my @v13             = splice(@_, 0, 3);
    my @m               = M_Transpose(splice(@_, 0, 9));
    my @v21             = splice(@m, 0, 3);
    my @v22             = splice(@m, 0, 3);
    my @v23             = splice(@m, 0, 3);

    return (
      V_Dot(@v11, @v21), V_Dot(@v11, @v22), V_Dot(@v11, @v23),
      V_Dot(@v12, @v21), V_Dot(@v12, @v22), V_Dot(@v12, @v23),
      V_Dot(@v13, @v21), V_Dot(@v13, @v22), V_Dot(@v13, @v23));
  }


  sub M_Unit { return (1,0,0, 0,1,0, 0,0,1); }

  sub PI { return 4*atan2(1,1); }

  sub M_Rotate
  {                                                     # vmd convention
    my $n               = shift(@_);
    my $alpha           = shift(@_)*PI()/180;
    my $cos             = cos($alpha);
    my $sin             = sin($alpha);

    $cos                = 0 if (abs($cos)<1e-16);
    $sin                = 0 if (abs($sin)<1e-16);
    return (1,0,0, 0,$cos,-$sin, 0,$sin,$cos) if ($n==0); # around x-axis
    return ($cos,0,$sin, 0,1,0, -$sin,0,$cos) if ($n==1); # around y-axis
    return ($cos,-$sin,0, $sin,$cos,0, 0,0,1) if ($n==2); # around z-axis
    return M_Unit();
  }


  sub MV_Dot
  {
    my @v11             = splice(@_, 0, 3);
    my @v12             = splice(@_, 0, 3);
    my @v13             = splice(@_, 0, 3);
    my @v2              = splice(@_, 0, 3);

    return (V_Dot(@v11, @v2), V_Dot(@v12, @v2), V_Dot(@v13, @v2));
  }


# CHARMM input

  sub PSFConnectivity
  {
    my $n               = PSFGoto(bonds);

    return if (scalar(@nconnect));
    printf("Info: creating connectivity\n") if ($info);
    for (my $i=0; $i<$n; ++$i)
    {
      my @bond          = PSFGet(2);
      $connect[$bond[0]][$nconnect[$bond[0]]++] = $bond[1];
      $connect[$bond[1]][$nconnect[$bond[1]]++] = $bond[0];
    }
  }


  sub PSFDihedrals                                      # hack to accomodate
  {                                                     # LAMMPS' way of calc
    $idihedral          = 0;                            # LJ 1-4 interactions
    return $ndihedral if (($dihedral_flag = $ndihedral ? 1 : 0));
    PSFConnectivity();
    printf("Info: creating dihedrals\n") if ($info);
    my $n               = scalar(@nconnect);
    my @bonded          = ();
    for (my $i=1; $i<=$n; ++$i)
    {
      $bonded[0]        = $i;
      for (my $i=0; $i<scalar($nconnect[$bonded[0]]); ++$i)
      {
        $bonded[1]      = $connect[$bonded[0]][$i];
        for (my $i=0; $i<scalar($nconnect[$bonded[1]]); ++$i)
        {
          next if (($bonded[2] = $connect[$bonded[1]][$i])==$bonded[0]);
          for (my $i=0; $i<scalar($nconnect[$bonded[2]]); ++$i)
          {
            next if (($bonded[3] = $connect[$bonded[2]][$i])==$bonded[1]);
            next if ($bonded[3]<$bonded[0]);
            $dihedral[$ndihedral++] = join(" ", @bonded);
          }
        }
      }
    }
    $dihedral_flag      = 1;
    return $ndihedral;
  }
        

  sub CreatePSFIndex                                    # make an index of id
  {                                                     # locations
    my @psf_ids         = ("!NATOM","!NBOND:","!NTHETA:","!NPHI:","!NIMPHI:");
    my @ids             = (atoms, bonds, angles, dihedrals, impropers);
    my $k               = 0;
    my %hash;

    printf("Info: creating PSF index\n") if ($info);
    open(PSF, "<$project.psf") if (fileno(PSF) eq "");
    foreach (@psf_ids) { $hash{$_} = shift(@ids); };
    while (<PSF>)
    {
      chop();
      my @tmp           = split(" ");
      my $n             = $hash{$tmp[1]};
      $PSFIndex{$n}     = tell(PSF)." ".$tmp[0] if ($n ne "");
    }
  }


  sub PSFGoto                                           # goto $ident in <PSF>
  {
    CreatePSFIndex() if (!scalar(%PSFIndex));
    my $id              = shift(@_);
    my @n               = split(" ", $PSFIndex{$id});

    @PSFBuffer          = ();
    # return PSFDihedrals() if ($id eq "dihedrals");
    if (!scalar(@n))
    {
      printf("Warning: PSF index for $id not found\n");
      seek(PSF, 0, SEEK_END);
      return -1;
    }
    seek(PSF, $n[0], SEEK_SET);
    return $n[1];
  }


  sub PSFGet
  {
    if ($dihedral_flag)
    {
      $dihedral_flag    = $idihedral+1<$ndihedral ? 1 : 0;
      return split(" ", $dihedral[$idihedral++]);
    }
    if (!scalar(@PSFBuffer))
    {
      my $line          = <PSF>;
      chop($line);
      @PSFBuffer        = split(" ", $line);
    }
    return splice(@PSFBuffer, 0, shift(@_));
  }


  sub PSFWrite
  {
    my $items           = shift(@_);
    my $n               = $items;

    if ($psf_ncols>7) { printf(PSF_CTRL "\n"); $psf_ncols = 0; }
    foreach(@_)
    {
      printf(PSF_CTRL " %7d", $_);
      ++$psf_ncols;
      if ((!--$n) && ($psf_ncols>7))
      {
        printf(PSF_CTRL "\n");
        $psf_ncols      = 0;
        $n              = $items;
      }
    }
  }


  sub CRDGoto
  {
    my $n;

    return if (shift(@_) ne "atoms");
    open(CRD, "<".($pdb ? $Pdb : $Crd)) if (fileno(CRD) eq "");
    seek(CRD, 0, SEEK_SET);
    return PSFGoto(atoms) if ($pdb);
    while (substr($n = <CRD>, 0, 1) eq "*") {}
    chop($n);
    return $n;
  }


  sub NextPDB2CRD
  {
    my @n = (6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,6,4,2,2);
    my @data = ();
    my $c = 0;
    my $line;

    while (substr($line = <CRD>, 0, 4) ne "ATOM") {};
    chop($line);
    foreach (@n) { push(@data, substr($line, ($c += $_)-$_, $_)); }
    return @data[1, 8, 5, 3, 11, 12, 13, 17, 8, 15];
  }


  sub Delete
  {
    my $item            = shift(@_);
    my $k               = 0;
    my @list;

    foreach (@_)
    {
      my @tmp           = split(" ");
      delete($tmp[$item]);
      $list[$k++]       = join(" ", @tmp);
    }
    return @list;
  }


  sub CreateID                                          # create id from list
  {
    my $n               = scalar(@_);
    my @list            = @_;
    my $id              = "";
    my $flag            = $list[0] gt $list[-1];
    my $j               = $n;
    my $tmp;

    return "" if (scalar(@list)<$n);
    $flag               = $list[1] gt $list[-2]
                          if ((scalar(@list)>3)&&($list[0] eq $list[-1]));
    for (my $i=0; $i<$n; ++$i)
    {
      $id               .= ($i ? " " : "").($tmp = $list[$flag ? --$j : $i]);
      $id               .= substr("    ", 0, 4-length($tmp));
    }
    return $id;
  }


  sub AtomTypes
  {
    my $n               = PSFGoto(atoms);
    my %list;

    return () if ($n<1);
    $atom_types[0]      = -1;
    for (my $i=0; $i<$n; ++$i)
    {
      my @tmp           = split(" ", <PSF>);
      $tmp[5]           = $symbols{$tmp[5]}
        if ((substr($tmp[5],0,1) lt '0')||(substr($tmp[5],0,1) gt '9'));
      push(@atom_types, $tmp[5]);
      ++$list{$tmp[5]};
    }
    if ($water_dens)
    {
      push(@atom_types, $symbols{HT}); ++$list{$symbols{HT}};
      push(@atom_types, $symbols{OT}); ++$list{$symbols{OT}};
    }
    if ($ions)
    {
      push(@atom_types, $symbols{CLA}); ++$list{$symbols{CLA}};
      push(@atom_types, $symbols{SOD}); ++$list{$symbols{SOD}};
    }
    return sort({$a<=>$b} keys(%list));
  }


  sub Markers
  {
    my %markers = (
      NONBONDED => '0',
      BONDS     => '1',
      ANGLES    => '2',
      DIHEDRALS => '3',
      IMPROPERS => '4',
      IMPROPER  => '4'
    );
    return %markers;
  }


  sub NonBond
  {
    my @cols            = @_;
    my $f               = (scalar(@cols)>3)&&(substr($cols[3],0,1) ne "!");
    my @tmp             = (-$cols[1], $cols[2],
                            $f ? -$cols[4]:-$cols[1], $f ? $cols[5]:$cols[2]);
    $tmp[1]             *= 2.0**(5/6);                  # adjust sigma
    $tmp[3]             *= 2.0**(5/6);                  # adjust sigma 1-4
    return join(" ", @tmp);
  }


  sub AtomParameters                                    # non-bonded parameters
  {
    my @types;
    my @list;
    my $k               = 0;
    my $read            = 0;
    my %markers         = Markers();

    foreach(@_) { $types{$ids{$_}} = $k++; }
    seek(PARAMETERS, 0, 0);
    while (<PARAMETERS>)
    {
      chop();
      my @cols          = split(" ");
      if ($read&&(scalar(@cols)>1)&&
        (substr($cols[0],0,1) ne "!")&&($cols[1] lt "A"))
      {
        my $k           = $types{shift(@cols)};
        $list[$k]       = NonBond(@cols) if ($k ne "");
      }
      if ($markers{$cols[0]} ne "") {
        $read           = ($markers{$cols[0]} eq "0") ? 1 : 0; }
    }
    $list[$types{HT}]   = NonBond(0, -0.046, 0.2245)
      if ($water_dens&&($list[$types{HT}] eq ""));
    $list[$types{OT}]   = NonBond(0, -0.152100, 1.768200)
      if ($water_dens&&($list[$types{OT}] eq ""));
    $list[$types{CLA}]  = NonBond(0, -0.150, 2.27)
      if ($ions&&($list[$types{CLA}] eq ""));
    $list[$types{SOD}] = NonBond(0, -0.0469, 1.36375)
      if ($ions&&($list[$types{SOD}] eq ""));
    return @list;
  }


  sub BondedTypes                                       # create bonded types
  {
    my $mode            = shift(@_);                    # operation mode
    my $items           = (2, 3, 4, 4)[$mode];          # items per entry
    my $id              = (bonds, angles, dihedrals, impropers)[$mode];
    my $n               = PSFGoto($id);
    my %list;

    for (my $i=0; $i<$n; ++$i)
    {
      my @tmp           = ();
      foreach (PSFGet($items)) { push(@tmp, $ids{$atom_types[$_]}); }
      ++$list{CreateID(@tmp)};
    }
    ++$list{CreateID(HT, OT)} if ($water_dens&&($mode==0));
    ++$list{CreateID(HT, OT, HT)} if ($water_dens&&($mode==1));
    @types              = sort(keys(%list));
  }


  sub Parameters                                        # parms from columns
  {
    my $items           = shift(@_);
    my @cols            = @_;
    my $parms           = "";

    for (my $i=$items; ($i<scalar(@cols))&&(substr($cols[$i],0,1)ne"!"); ++$i)
    {
      $parms            = $parms.($i>$items ? " " : "").$cols[$i];
    }
    return $parms;
  }


  sub BondedParameters                                  # distil parms from
  {                                                     # <PARAMETERS>
    my $mode            = shift(@_);                    # bonded mode
    return if (($mode>3)||($mode<0));

    my $items           = (2, 3, 4, 4)[$mode];          # items per entry
    my $name            = ("bond", "angle", "dihedral", "improper")[$mode];
    my $read            = 0;
    my $k               = 0;
    my %markers         = Markers();
    my @set;
    my @tmp;
    my $f;
    my %list;
    my %link;

    @parms              = ();
    foreach(@types) { $link{$_} = $k++; }
    seek(PARAMETERS, 0, 0);
    while (<PARAMETERS>)
    {
      chomp();
      my @cols          = split(" ");
      if ($read&&(scalar(@cols)>$items)&&($cols[$items] lt "A"))
      {
        if (($items==4)&&(($f = ($cols[1] eq "X")&&($cols[2] eq "X"))||
          (($cols[0] eq "X")&&($cols[3] eq "X"))))      # wildcards
        {
          my $id        = CreateID(($cols[1-$f], $cols[2+$f]));
          for ($k=0; $k<scalar(@types); ++$k)
          {
            if (!$set[$k])
            {
              my @tmp   = split(" ", $types[$k]);
              if (CreateID($tmp[1-$f], $tmp[2+$f]) eq $id)
              {
                if ($mode==2)
                {
                  if ($parms[$k] eq "") {
                    $parms[$k] = Parameters($items,@cols)." 1"; }
                  else {
                    $parms[$k] .= ":".Parameters($items,@cols)." 0"; }
                }
                else {
                  $parms[$k] .= Parameters($items,@cols); }
              }
            }
          }
        }
        else                                            # regular
        {
          for (my $i=0; $i<$items; ++$i) { $tmp[$i] = $cols[$i]; };
          $k                    = $link{CreateID(@tmp)};
          if ($k ne "")
          {
            $parms[$k]  = "" if (!$set[$k]);
            $parms[$k]  .= ($set[$k]++ ? ":" : "").Parameters($items,@cols);
            $parms[$k]  .= ($set[$k]-1 ? " 0" : " 1") if ($mode==2);
          }
        }
      }
      if ($markers{$cols[0]}) {
        $read           = ($markers{$cols[0]} eq $mode+1) ? 1 : 0; }
    }
    if ($water_dens)
    {
      $parms[$link{CreateID(HT, OT)}] = "450 0.9572" if ($mode==0);
      $parms[$link{CreateID(HT, OT, HT)}] = "55 104.52" if ($mode==1);
    }
    for (my $i=0; $i<scalar(@types); ++$i)
    {
      printf("Warning: %s parameter %4d for [%s] was not found\n",
        $name, $i+1, $types[$i]) if ($parms[$i] eq "");
    }
  }


  sub SetScreeningFactor                                # set screening factor
  {
    my $id              = shift(@_);
    my $value           = shift(@_);
    my $new             = "";

    foreach (split(":", $parms[$id]))
    {
      my @tmp           = split(" ");
      $tmp[-1]          = $value if ($tmp[-1]);
      $new              .= ":" if ($new ne "");
      $new              .= join(" ", @tmp);
    }
    $parms[$id]         = $new;
  }


  sub CorrectDihedralParameters
  {
    my $n               = PSFGoto(dihedrals);
    my %hash;
    my $hash_id;
    my $id1;
    my $id2;
    my $first;
    my $last;

    for (my $i=0; $i<$n; ++$i)
    {
      my @bonded        = PSFGet(4);
      my @tmp           = ();
      foreach (@bonded) { push(@tmp, $ids{$atom_types[$_]}); }
      $id1              = $link{CreateID(@tmp)}-1;
      $first            = $bonded[0];
      $last             = $bonded[3];
      if ($first>$last) { my $tmp = $first; $first = $last; $last = $tmp; }
      if (($id2 = $hash{$hash_id = $first." ".$last}) eq "")
      {
        $hash{$hash_id} = $id1;                         # add id to hash
      }
      else
      {
        SetScreeningFactor($id1, 0.5);                  # 6-ring: shared 1-4
        SetScreeningFactor($id2, 0.5);
      }
    }
    $n                  = PSFGoto(angles);
    for (my $i=0; $i<$n; ++$i)
    {
      my @bonded        = PSFGet(3);
      $first            = $bonded[0];
      $last             = $bonded[2];
      if ($first>$last) { my $tmp = $first; $first = $last; $last = $tmp; }
      if (($id1 = $hash{$first." ".$last}) ne "")
      {
        SetScreeningFactor($id1, 0);                    # 5-ring: no 1-4
      }
    }
  }


  sub AddMass
  {
    my $symbol          = shift(@_);
    my $mass            = shift(@_);

    return if ($symbols{$symbol} ne "");
    $ids{++$max_id}     = $symbol;
    $masses{$max_id}    = $mass;
    $symbols{$symbol}   = $max_id;
  }


  sub ReadTopology                                      # read topology links
  {
    my $id              = shift(@_);
    my $item            = shift(@_);
    my $read            = 0;
    my @tmp;

    open(TOPOLOGY, "<top_$forcefield.rtf");
    $max_id             = 0;
    while (<TOPOLOGY>)
    {
      chop();                                           # delete CR at end
      my @tmp           = split(" ");
      $read             = 1 if ($tmp[0] eq "MASS");
      if ($read&&($tmp[0] eq "MASS"))
      {
        $symbols{$tmp[2]}       = $tmp[1];
        $ids{$tmp[1]}           = $tmp[2];
        $masses{$tmp[1]}        = $tmp[3];
        $max_id                 = $tmp[1] if ($max_id<$tmp[1]);
      } elsif ($read&&($tmp[0] eq "ATOM")) {
        # quit reading when hitting the "ATOM" section
        last;
      }
    }
    AddMass(HT, 1.00800);
    AddMass(OT, 15.99940);
    AddMass(CLA, 35.450000);
    AddMass(SOD, 22.989770);
    close(TOPOLOGY);
  }


  sub CrossLink                                         # symbolic cross-links
  {
    my @list            = @_;
    my $n               = scalar(@list);
    my %hash;

    for (my $i=0; $i<$n; ++$i) { $hash{$list[$i]} = $i+1; }
    return %hash;
  }


  sub CharacterizeBox
  {
    my $flag            = 1;
    my @x               = (-$L[0]/2, $L[0]/2);
    my @y               = (-$L[1]/2, $L[1]/2);
    my @z               = (-$L[2]/2, $L[2]/2);
    my $n               = CRDGoto(atoms);
    my $extremes        = !($L[0] && $L[1] && $L[2]);

    @Center             = (0, 0, 0);
    return if (!$n);
    for (my $i=0; $i<$n; ++$i)
    {
      my @tmp           = $pdb ? NextPDB2CRD() : split(" ", <CRD>);
      my @p             = @tmp[-6, -5, -4];
      @p                = MV_Dot(@R, @p);
      $x[0]             = $p[0] if ($flag||($p[0]<$x[0]));
      $x[1]             = $p[0] if ($flag||($p[0]>$x[1]));
      $y[0]             = $p[1] if ($flag||($p[1]<$y[0]));
      $y[1]             = $p[1] if ($flag||($p[1]>$y[1]));
      $z[0]             = $p[2] if ($flag||($p[2]<$z[0]));
      $z[1]             = $p[2] if ($flag||($p[2]>$z[1]));
      $flag             = 0 if ($flag);
    }
    $L[0]               = $x[1]-$x[0] if (!$L[0]);
    $L[1]               = $y[1]-$y[0] if (!$L[1]);
    $L[2]               = $z[1]-$z[0] if (!$L[2]);
    $L[0]               += $border;
    $L[1]               += $border;
    $L[2]               += $border;
    @Center             = (($x[1]+$x[0])/2, ($y[1]+$y[0])/2, ($z[1]+$z[0])/2);
    printf("Info: recentering atoms\n") if ($info&&$center);
  }


  sub SetupWater
  {
    return if (!$water_dens);

    my $dens            = 1000*$water_dens;             # kg/m^3
    my $m               = 0.018;                        # kg/mol
    my $loh             = 0.9572;                       # l[O-H] in [A]
    my $s_OT            = 1.7682;                       # CHARMM sigma [A]
    my $ahoh            = (180-104.52)/360*PI();
    my @p               = ($loh*cos($ahoh), $loh*sin($ahoh), 0);

    printf("Info: creating fcc water\n") if ($info);
    $n_water            = 4;                            # molecules/cell
    $nav                = 6.022e23;                     # 1/mol
    $v_water            = $m/$nav/$dens*1e30;           # water volume [A^3]
    $r_water            = $s_OT*2**(-1/6);              # sigma_OT in [A]
    @p_water            = (0,0,0, @p, -$p[0],$p[1],0);
    $v_fcc              = $n_water*$v_water;            # cell volume
    $l_fcc              = $v_fcc**(1/3);                # cell length
    @p_fcc              = (0.00,0.00,0.00, 0.50,0.50,0.00,
                           0.50,0.00,0.50, 0.00,0.50,0.50);
    @n_fcc              = ();
    for (my $i=0; $i<scalar(@L); ++$i)
    {
      my $n             = $L[$i]/$l_fcc;                # calculate n_fcc
      $n                = int($n-int($n) ? $n+1 : $n);  # ceil($n)
      $L[$i]            = $n*$l_fcc;                    # adjust box length
      printf("Info: changed l%s to %g A\n", ("x","y","z")[$i], $L[$i])
        if ($info);
      push(@n_fcc, $n);
    }
    foreach (@p_fcc) { $_ = ($_+0.25)*$l_fcc; }         # p_fcc in [A]
    for (my $x=0; $x<$n_fcc[0]; ++$x) {                 # initialize flags
      for (my $y=0; $y<$n_fcc[1]; ++$y) {
        for (my $z=0; $z<$n_fcc[2]; ++$z) {
          $flags_fcc[$x][$y][$z] = 15; } } }            # turn on all fcc sites
  }


  sub floor
  {
    my $x               = shift(@_);

    return $x>0 ? int($x) : int($x)-1;
  }


  sub Periodic
  {
    my @p               = splice(@_, 0, 3);

    return (
      $p[0]-floor($p[0]/$L[0]+0.5)*$L[0],
      $p[1]-floor($p[1]/$L[1]+0.5)*$L[1],
      $p[2]-floor($p[2]/$L[2]+0.5)*$L[2]);
  }


  sub EraseWater
  {
    my $r               = shift(@_)/2;
    my @p               = splice(@_, 0, 3);
    @p                  = V_Subtr(@p, @Center) if (!$center);
    my @edges           = (
      $p[0]-$r,$p[1]-$r,$p[2]-$r, $p[0]-$r,$p[1]-$r,$p[2]+$r,
      $p[0]-$r,$p[1]+$r,$p[2]-$r, $p[0]-$r,$p[1]+$r,$p[2]+$r,
      $p[0]+$r,$p[1]-$r,$p[2]-$r, $p[0]+$r,$p[1]-$r,$p[2]+$r,
      $p[0]+$r,$p[1]+$r,$p[2]-$r, $p[0]+$r,$p[1]+$r,$p[2]+$r);
    my %list;
    my @n;

    my $d2              = ($r_water+$r)**2;
    my @l               = ($L[0]/2, $L[1]/2, $L[2]/2);
    for (my $i=0; $i<scalar(@edges); $i+=3)             # determine candidates
    {
      my @q             = Periodic(@edges[$i, $i+1, $i+2]);
      my @n             = (int(($q[0]+$l[0])/$l_fcc),int(($q[1]+$l[1])/$l_fcc),
                           int(($q[2]+$l[2])/$l_fcc));
      ++$list{join(" ", @n)};
    }
    foreach (sort(keys(%list)))                         # check overlap
    {
      my @n             = split(" ");
      my @corner        = ($n[0]*$l_fcc-$l[0]+$p_water[0],
                           $n[1]*$l_fcc-$l[1]+$p_water[1],
                           $n[2]*$l_fcc-$l[2]+$p_water[2]);
      my $bit           = 1;
      my $flags         = 0;
      for (my $i=0; $i<scalar(@p_fcc); $i+=3)
      {
        my @q           = V_Add(@corner, @p_fcc[$i,$i+1,$i+2]);
        my @dp          = Periodic(V_Subtr(@q, @p));
        $flags          |= $bit if (V_Dot(@dp, @dp)>$d2); # turn on fcc
        $bit            *= 2;
      }
      $flags_fcc[$n[0]][$n[1]][$n[2]] &= $flags;        # set flags
    }
  }


  sub CountFCC
   {
    my $n               = 0;

    return $n_fccs = 0 if (!$water_dens);
    for (my $x=0; $x<$n_fcc[0]; ++$x) {                 # count water
      for (my $y=0; $y<$n_fcc[1]; ++$y) {
        for (my $z=0; $z<$n_fcc[2]; ++$z) {
          my $bit       = 1;
          my $flags     = $flags_fcc[$x][$y][$z];
          for (my $i=0; $i<$n_water; ++$i) {
            ++$n if ($flags & $bit);
            $bit        *= 2; } } } }
    return ($n_fccs = $n);
  }


  sub AddIons
  {
    my $n               = ($n_waters = CountFCC())-int(abs($net_charge));

    return if (!$ions);
    printf("Warning: charge not neutralized: too little water\n") if ($n<0);
    return if ($n<0);
    printf(
      "Warning: charge not neutralized: net charge (%g) is not an integer\n",
      $net_charge) if ($net_charge!=int($net_charge));
    my $n_na            = $net_charge<0 ? int(abs($net_charge)) : 0;
    my $n_cl            = $net_charge>0 ? int(abs($net_charge)) : 0;
    my $n_mol           = int($ion_molar*$n*$v_water*1e-27*$nav+0.5);
    my $n_atoms         = ($n_na += $n_mol)+($n_cl += $n_mol);
    $n_waters           -= $n_atoms;
    printf(
      "Info: adding ions: [NaCl] = %g mol/l (%d Na+, %d Cl-)\n",
      $n_mol/$n/$v_water/$nav/1e-27, $n_na, $n_cl) if ($info);
    $n                  += int(abs($net_charge));
    my $salt            = 2**$n_water;
    srand(time());                                      # seed random number
    for (my $x=0; $x<$n_fcc[0]; ++$x)                   # replace water by ions
    {
      for (my $y=0; $y<$n_fcc[1]; ++$y)
      {
        for (my $z=0; $z<$n_fcc[2]; ++$z)
        {
          my $bit       = 1;
          my $flags     = $flags_fcc[$x][$y][$z];
          for (my $i=0; $i<$n_water; ++$i)
          {
            if ($flags & $bit)
            {
              my $prob  = $n_atoms/$n;
              --$n;
              if (rand()<$prob)
              {
                my $na  = rand()<$n_na/$n_atoms ? 1 : 0;
                --$n_atoms;
                if ($na) { --$n_na; } else { --$n_cl; }
                $flags  |= $salt*(1+$salt*$na)*$bit;    # set type of ion
              }
            };
            $bit        *= 2;
          }
          $flags_fcc[$x][$y][$z] = $flags;
        }
      }
    }
  }


# LAMMPS output

  sub WriteLAMMPSHeader                                 # print lammps header
  {
    printf(LAMMPS "LAMMPS data file. %sCreated by $program v$version on %s\n",
           ($add ? "CGCMM Style. atom_style full. " : ""),`date`);
    printf(LAMMPS "%12d  atoms\n", $natoms);
    printf(LAMMPS "%12d  bonds\n", $nbonds);
    printf(LAMMPS "%12d  angles\n", $nangles);
    printf(LAMMPS "%12d  dihedrals\n", $ndihedrals);
    printf(LAMMPS "%12d  impropers\n\n", $nimpropers);
    printf(LAMMPS "%12d  atom types\n", $natom_types);
    printf(LAMMPS "%12d  bond types\n", $nbond_types);
    printf(LAMMPS "%12d  angle types\n", $nangle_types);
    printf(LAMMPS "%12d  dihedral types\n", $ndihedral_types);
    printf(LAMMPS "%12d  improper types\n\n", $nimproper_types);
  }


  sub WriteControlHeader
  {
    printf(PDB_CTRL "REMARK  \n");
    printf(PDB_CTRL "REMARK  CONTROL PDB %s_ctrl.pdb\n", $project);
    printf(PDB_CTRL "REMARK  CREATED BY %s v%s ON %s",
      $program, $version, `date`);
    printf(PDB_CTRL "REMARK  \n");

    printf(PSF_CTRL "PSF\n");
    printf(PSF_CTRL "\n");
    printf(PSF_CTRL "%8d !NTITLE\n", 2);
    printf(PSF_CTRL " REMARKS CONTROL PSF %s_ctrl.psf\n", $project);
    printf(PSF_CTRL " REMARKS CREATED BY %s v%s ON %s",
      $program, $version, `date`);
    printf(PSF_CTRL "\n");
  }


  sub WriteBoxSize                              # print box limits
  {
    my @lo              = V_Mult(@L[0,1,2], -1/2);
    my @hi              = V_Mult(@L[0,1,2], 1/2);

    @lo                 = V_Add(@lo, @Center) if (!$center);
    @hi                 = V_Add(@hi, @Center) if (!$center);
    printf(LAMMPS "%12.8g %12.8g xlo xhi\n", $lo[0], $hi[0]);
    printf(LAMMPS "%12.8g %12.8g ylo yhi\n", $lo[1], $hi[1]);
    printf(LAMMPS "%12.8g %12.8g zlo zhi\n\n", $lo[2], $hi[2]);
  }


  sub WriteMasses                                       # print mass list
  {
    my $k               = 0;

    printf(LAMMPS "Masses\n\n");
    foreach (@types)
    {
      printf(LAMMPS "%8d %10.7g%s\n",
        ++$k, $masses{$_}, $add ? "  # ".$ids{$_} : "");
    }
    printf(LAMMPS "\n");
  }


  sub WriteFCCAtoms
  {
    my $k               = shift(@_);
    my $res             = shift(@_);

    return $k if (!$water_dens);
    $k_fcc              = $k+1;
    my @id              = ($symbols{OT}, $symbols{HT}, $symbols{HT},
                           $symbols{SOD}, $symbols{CLA});
    my @par             = ();
    my @charge          = (-0.834, 0.417, 0.417, 1, -1);
    my $salt            = 2**$n_water;
    my @l               = ($L[0]/2, $L[1]/2, $L[2]/2);
    my $iwater          = 0;
    my $isalt           = 0;
    foreach(@id) { push(@par, $link{$_}); }
    for (my $x=0; $x<$n_fcc[0]; ++$x)
    {
      for (my $y=0; $y<$n_fcc[1]; ++$y)
      {
        for (my $z=0; $z<$n_fcc[2]; ++$z)
        {
          my @corner    = ($x*$l_fcc-$l[0], $y*$l_fcc-$l[1], $z*$l_fcc-$l[2]);
          my $flags     = $flags_fcc[$x][$y][$z];
          my $bit       = 1;
          for (my $i=0; $i<scalar(@p_fcc); $i+=3)
          {
            my $pair    = $bit;
            if ($flags & $pair)
            {
              my @p     = V_Add(@corner, @p_fcc[$i,$i+1,$i+2]);
              my $j     = 0;                            # print water
              my $n     = scalar(@p_water);
              ++$res;
              if ($flags & ($pair *= $salt))            # print salt ion
              {                                         # sodium if highest
                $j      = $flags & ($pair*$salt) ? 3 : 4;
                $n      = 1;
                $counter = ++$isalt;
              }
              else { $counter = ++$iwater; }
              for (my $i=0; $i<$n; $i+=3)
              {
                my @xyz = V_Add(@p, @p_water[$i,$i+1,$i+2]);
                @xyz    = V_Add(@xyz, @Center) if (!$center);
                printf(LAMMPS "%8d %7d %5d %9.6g %16.12g %16.12g %16.12g%s\n",
                  ++$k, $res, $par[$j], $charge[$j], $xyz[0], $xyz[1],
                  $xyz[2], $add ? "  # ".$types[$par[$j]-1] : "");
                printf(PDB_CTRL "ATOM %6.6s %-4.4s %-3.3s %5.5s %3.3s ".
                  "%7.7s %7.7s %7.7s %5.5s %5.5s %4.4s %s\n", $k,
                  $types[$par[$j]-1], $n-1 ? "HOH" : "ION", $res, "",
                  $xyz[0], $xyz[1], $xyz[2], "1.00", "0.00", "",
                  $n-1 ? "WATR" : "SALT") if ($pdb_ctrl);
                printf(PSF_CTRL "%8d %4.4s %-4.4s %-4.4s %-4.4s %4.4s ".
                  "%16.8e %7.7s %9.9s 0\n", $k, $n-1 ? "WATR" : "SALT",
                  $counter, $n-1 ? "HOH" : "ION", $types[$par[$j]-1], $id[$j],
                  $charge[$j], $masses{$id[$j]}, "") if ($pdb_ctrl);
                ++$j;
              }
            }
            $bit        *= 2;
          }
        }
      }
    }
    return $k;
  }


  sub WritePSFAtoms()
  {
    my $n               = PSFGoto(atoms);
    my @res             = (0, 0);

    printf(PSF_CTRL "%8d !NATOM\n", $n+2*$n_waters+$n_fccs);
    while (<PSF>)
    {
      last if (!$n--);
      my @psf           = split(" ");
      if ($res[1]!=$psf[2]) { ++$res[0]; $res[1] = $psf[2]; }
      printf(PSF_CTRL "%8d %4.4s %-4.4s %-4.4s %-4.4s %-4.4s ".
        "%16.8e %7.7s %9.9s %s\n", $psf[0], $psf[1], $res[0],
        $psf[3], $psf[4], $psf[5], $psf[6], $psf[7], "", $psf[8]);
    }
  }


  sub WriteAtoms                                        # print positions etc.
  {
    my $n               = PSFGoto(atoms);
    my $k               = 0;
    my @res             = (0, 0);

    CRDGoto(atoms);
    $net_charge         = 0;
    printf(LAMMPS "Atoms%s\n\n",($add ? " # full" : "")) if ($n>0);
    for (my $i=0; $i<$n; ++$i)
    {
      my @crd           = $pdb ? NextPDB2CRD() : split(" ", <CRD>);
      my @psf           = split(" ", <PSF>);
      my @xyz           = MV_Dot(@R, @crd[-6, -5, -4]);
      @xyz              = V_Subtr(@xyz, @Center) if ($center);
      if ($crd[-2]!=$res[1]) { ++$res[0]; $res[1] = $crd[-2]; }
      printf(LAMMPS "%8d %7d %5d %9.6g %16.12g %16.12g %16.12g%s\n", ++$k,
        $res[0], $link{$atom_types[$k]}, $psf[6], $xyz[0], $xyz[1], $xyz[2],
        $add ? "  # ".$types[$link{$atom_types[$k]}-1] : "");
      printf(PDB_CTRL "ATOM %6.6s %-4.4s %-4.4s %4.4s %3.3s ".
        "%7.7s %7.7s %7.7s %5.5s %5.5s %4.4s %s\n", $k,
        $crd[-7], $crd[-8], $res[0], "", $xyz[0], $xyz[1], $xyz[2],
        "1.00", $crd[-1], "", $crd[-3]) if ($pdb_ctrl);
      next if (!$water_dens);                           # is water added?
      $net_charge       += $psf[6];
      my @c             = split(" ", $parms[$link{$atom_types[$k]}-1]);
      EraseWater($c[1], @xyz);
    }
    $net_charge         = int($net_charge*1e5+($net_charge>0?0.5:-0.5))/1e5;
    AddIons() if ($water_dens);
    WritePSFAtoms() if ($pdb_ctrl);
    $k                  = WriteFCCAtoms($k, $res[0]+$res[1]);
    printf(PDB_CTRL "END\n") if ($pdb_ctrl);
    printf(LAMMPS "\n");
    return $k;
  }


  sub WriteParameters                           # print parameters
  {
    my $mode            = shift(@_)+1;
    my $header          = ("Pair","Bond","Angle","Dihedral","Improper")[$mode];
    my $hint            = ("# lj/charmm/coul/long", "# harmonic", "# charmm", "# charmm", "# harmonic")[$mode];
    my $n               = (4, 2, 4, 4, 2)[$mode];
    my $k               = 0;

    printf("Info: converting ".lc($mode ? $header : "Atom")."s\n") if ($info);
    if ($mode--)
    {
      BondedTypes($mode);
      BondedParameters($mode);
      %link             = CrossLink(@types);
      CorrectDihedralParameters() if ($mode==2);
      @parms            = Delete(1, @parms) if ($mode==3);
    }
    return 0 if (!scalar(@parms));
    printf(LAMMPS "%s Coeffs %s\n\n", $header, ($add ? $hint : ""));
    for (my $i=0; $i<scalar(@parms); ++$i)
    {
      if ($parms[$i] ne "")
      {
        foreach (split(":", $parms[$i]))
        {
          my @tmp       = split(" ");
          printf(LAMMPS "%8d", ++$k);
          for (my $j=0; $j<$n; ++$j) {
            printf(LAMMPS " %16.12g", $j<scalar(@tmp) ? $tmp[$j] : 0); }
          printf(LAMMPS "%s\n", $add ? "  # ".$types[$i] : "");
        }
      } else { ++$k; }
    }
    printf(LAMMPS "\n");
    return $k;
  }


  sub WriteFCCBonded
  {
    my $mode            = shift(@_);
    my $k               = shift(@_);
    my $atom            = $k_fcc;

    return $k if (($mode>1)||!$water_dens);
    my $type            = $mode ? CreateID(HT, OT, HT) : CreateID(HT, OT);
    my $id              = $link{$type};
    my $salt            = 2**$n_water;
    for (my $x=0; $x<$n_fcc[0]; ++$x)
    {
      for (my $y=0; $y<$n_fcc[1]; ++$y)
      {
        for (my $z=0; $z<$n_fcc[2]; ++$z)
        {
          my @corner    = ($x*$l_fcc-$L[0]/2, $y*$l_fcc-$L[1]/2,
                           $z*$l_fcc-$L[2]/2);
          my $flags     = $flags_fcc[$x][$y][$z];
          my $bit       = 1;
          for (my $i=0; $i<scalar(@p_fcc); $i+=3)
          {
            if ($flags&$bit)
            {
              if ($flags&($bit*$salt)) { ++$atom; }
              else
              {
                printf(LAMMPS "%8d %7d %7d %7d%s\n", ++$k, $id, $atom,
                  $atom+1, $add ? "  # ".$type : "") if (!$mode);
                printf(LAMMPS "%8d %7d %7d %7d%s\n", ++$k, $id, $atom,
                  $atom+2, $add ? "  # ".$type : "") if (!$mode);
                printf(LAMMPS "%8d %7d %7d %7d %7d%s\n", ++$k, $id, $atom+1,
                  $atom, $atom+2, $add ? "  # ".$type : "") if ($mode);
                if ($pdb_ctrl)
                {
                  PSFWrite(2, $atom, $atom+1, $atom, $atom+2) if (!$mode);
                  PSFWrite(3, $atom+1, $atom, $atom+2) if ($mode);
                }
                $atom   += 3;
              }
            }
            $bit        *= 2;
          }
        }
      }
    }
    return $k;
  }


  sub WriteBonded                                       # print bonded list
  {
    my $mode            = shift(@_);
    my $psf_id          = ("!NBOND:", "!NTHETA:", "!NPHI:", "!NIMPHI:")[$mode];
    my $title           = ("bonds", "angles", "dihedrals", "impropers")[$mode];
    my $items           = (2, 3, 4, 4)[$mode];
    my $n               = PSFGoto($title);
    my $k               = 0;
    my @delta;
    my @tmp;

    return 0 if ($n<1);
    printf(LAMMPS "%s\n\n", ucfirst($title));
    printf(PSF_CTRL "\n%8d %s %s\n", $n+($mode ? ($mode==1 ? $n_waters : 0)
        : 2*$n_waters), $psf_id, $title) if ($pdb_ctrl);
    $psf_ncols          = 0 if ($pdb_ctrl);
    foreach (@parms)
    {
      push(@delta, $k);
      $k                += scalar(split(":"))-1 if ($_ ne "");
    }
    $k                  = 0;
    for (my $i=0; $i<$n; ++$i)
    {
      my @bonded        = PSFGet($items);
      my @tmp           = ();
      foreach (@bonded) { push(@tmp, $ids{$atom_types[$_]}); }
      my $id            = $link{CreateID(@tmp)}-1;
      my $m             = 0;
      if ($parms[$id] ne "")
      {
        foreach (split(":", $parms[$id]))
        {
          ++$m;
          my @const     = split(" ");
          next if (($const[0]==0)&&($mode==2 ? $const[-1]==0 : 1));
          printf(LAMMPS "%8d %7d", ++$k, $id+$delta[$id]+$m);
          foreach (@bonded) { printf(LAMMPS " %7d", $_); }
          printf(LAMMPS "%s\n", $add ? "  # ".CreateID(@tmp) : "");
        }
      }
      else
      {
        printf(LAMMPS "%8d %7d", ++$k, $id+$delta[$id]+$m);
        foreach (@bonded) { printf(LAMMPS " %7d", $_); }
        printf(LAMMPS "%s\n", $add ? "  # ".CreateID(@tmp) : "");
      }
      PSFWrite($items, @bonded) if ($pdb_ctrl);
    }
    $k                  = WriteFCCBonded($mode, $k);
    printf(PSF_CTRL "\n") if ($pdb_ctrl && $psf_ncols);
    printf(LAMMPS "\n");
    return $k;
  }


  sub CreateCorrectedPairCoefficients
  {
    my $read            = 0;
    my $k               = 0;
    my %id;
    my %type;

    $coefficients       = "";
    foreach (@types) { $id{$ids{$_}} = $_; $type{$_} = ++$k; }
    seek(PARAMETERS, 0, 0);
    while (<PARAMETERS>)
    {
      chop();
      my @cols          = split(" ");
      if ($read&&(scalar(@cols)>3)&&
        (substr($cols[0],0,1) ne "!")&&($cols[2] lt 'A'))
      {
        my $id1         = $id{$cols[0]};
        my $id2         = $id{$cols[1]};
        if (($id1 ne "")&&($id2 ne ""))
        {
          my @c         = (abs($cols[2]), $cols[3]*2.0**(-1/6));
          if ($type{$id2}<$type{$id1})
          {
            my $tmp = $id1; $id1 = $id2; $id2 = $tmp;
          }
          $coefficients .= ":" if ($coefficients ne "");
          $coefficients .= $type{$id1}." ".$type{$id2}." ";
          $coefficients .= $c[0]." ".$c[1]." ".$c[0]." ".$c[1];
        }
      }
      $read             = 1 if ($cols[0] eq "NBFIX");
      last if ($read&&!scalar(@cols));
    }
  }


  sub WriteData
  {
    open(LAMMPS, ">$project.in");                       # use .in for temporary
    open(PDB_CTRL, ">".$project."_ctrl.pdb") if ($pdb_ctrl);
    open(PSF_CTRL, ">".$project."_ctrl.psf") if ($pdb_ctrl);
    WriteControlHeader() if ($pdb_ctrl);
    ReadTopology();
    CharacterizeBox();
    SetupWater() if ($water_dens);
    WriteBoxSize();                             # body storage
    @types              = AtomTypes();                  # atoms
    @parms              = AtomParameters(@types);
    WriteMasses();
    %link               = CrossLink(@types);
    CreateCorrectedPairCoefficients();
    for (my $i=0; $i<scalar(@types); ++$i) { $types[$i] = $ids{$types[$i]}; }
    $natom_types        = WriteParameters(-1);  # pairs
    if ($#types+1 > $natom_types) {
      print "Warning: $#types atom types present, but only $natom_types pair coeffs found\n";
      # reset to what is found while determining the number of atom types.
      $natom_types = $#types+1;
    }
    $natoms             = WriteAtoms();
    $nbond_types        = WriteParameters(0);   # bonds
    $nbonds             = WriteBonded(0);
    $nangle_types       = WriteParameters(1);   # angles
    $nangles            = WriteBonded(1);
    $shake              = $link{CreateID(("HT", "OT", "HT"))};
    $ndihedral_types    = WriteParameters(2);   # dihedrals
    $ndihedrals         = WriteBonded(2);
    $nimproper_types    = WriteParameters(3);   # impropers
    $nimpropers         = WriteBonded(3);
    close(LAMMPS);                                      # close temp file
    open(LAMMPS, ">$project.data");                     # open data file
    WriteLAMMPSHeader();                                # header
    open(TMP, "<$project.in");                          # open temp file
    while (<TMP>) { printf(LAMMPS "%s", $_); }          # spool body
    close(TMP);                                         # close temp file
    if ($pdb_ctrl)
    {
      #while (<PSF>) { printf(PSF_CTRL "%s", $_); }
      close(PSF_CTRL); close(PDB_CTRL);
    }
    close(LAMMPS);                                      # close data file
  }


  sub WriteLAMMPSInput
  {
    open(LAMMPS, ">$project.in");                       # input file
    printf(LAMMPS "# Created by $program v$version on %s", `date`);
    printf(LAMMPS "# Command: %s\n\n", $cmdline);
    printf(LAMMPS "units           real\n");            # general
    printf(LAMMPS "neigh_modify    delay 2 every 1\n\n");
    printf(LAMMPS "atom_style      full\n");            # styles
    printf(LAMMPS "bond_style      harmonic\n") if ($nbond_types);
    printf(LAMMPS "angle_style     charmm\n") if ($nangle_types);
    printf(LAMMPS "dihedral_style  charmm\n") if ($ndihedral_types);
    printf(LAMMPS "improper_style  harmonic\n\n") if ($nimproper_types);
    printf(LAMMPS "pair_style      lj/charmm/coul/long 8 12\n");
    printf(LAMMPS "pair_modify     mix arithmetic\n");
    printf(LAMMPS "kspace_style    pppm 1e-6\n\n");

    if ($cmap) {
      printf(LAMMPS "# Modify following line to point to the desired CMAP file\n");
      printf(LAMMPS "fix             cmap all cmap charmm$cmap.cmap\n");
      printf(LAMMPS "fix_modify      cmap energy yes\n");
      printf(LAMMPS "read_data       $project.data fix cmap crossterm CMAP\n\n");
    }else{
      printf(LAMMPS "read_data       $project.data\n\n");       # read data
    }

    if ($coefficients ne "")                            # corrected coeffs
    {
      foreach (split(":", $coefficients))
      {
        printf(LAMMPS "pair_coeff      %s\n", $_);
      }
      printf(LAMMPS "\n");
    }
    printf(LAMMPS "special_bonds   charmm\n");          # invoke charmm
    printf(LAMMPS "thermo          10\n");              # set thermo style
    printf(LAMMPS "thermo_style    multi\n");
    printf(LAMMPS "timestep        1.0\n\n");           # 1.0 ps time step
    printf(LAMMPS "minimize 0.0 0.0 50 200\n\n");       # take of the edge
    printf(LAMMPS "reset_timestep  0\n");
    printf(LAMMPS "fix             1 all nve\n");
    printf(LAMMPS "fix             2 all shake 1e-6 500 0 m 1.0\n")
      if ($shake eq "");                                # shake all H-bonds
    printf(LAMMPS "fix             2 all shake 1e-6 500 0 m 1.0 a %s\n",$shake)
      if ($shake ne "");                                # add water if present
    printf(LAMMPS "velocity        all create 0.0 12345678 dist uniform\n\n");
    printf(LAMMPS "restart         500 $project.restart1 $project.restart2\n");
    printf(LAMMPS "dump            1 all atom 100 $project.dump\n");
    printf(LAMMPS "dump_modify     1 image yes scale yes\n\n");
    printf(LAMMPS "thermo          100\n");		# set thermo style
    printf(LAMMPS "run             1000\n");		# run for 1000 time steps
    close(LAMMPS);
  }

  # ----------------------- DESCRIPTION: sub CharmmCmap ------------------------ #
  # This subroutine add a new section "CMAP" to the LAMMPS data file, which      #
  # a part of the implementation of "CHARMM CMAP" (see references) in LAMMPS.    #
  # The section "CMAP" contains a list of dihedral ID pairs from adjecent        #
  # peptide backtone dihedrals whose dihedral angles are corrresponding to PHI   #
  # and PSI. (PHI: C--N--C_aphla_C and PSI: N--C_alpha--C--N)                    #
  #                                                                              #
  # Initiated by:   Xiaohu Hu (hux2@ornl.gov)                                    #
  # May 2009                                                                     #
  #                                                                              #
  # Finalized Oct 2016 by: Robert Latour (latourr@clemson.edu),                  #
  #   David Hyde-Volpe, and Tigran Abramyan, Clemson University,                 #
  #   and Chris Lorenz (chris.lorenz@kcl.ac.uk                                   #
  #                                                                              #
  # References:                                                                  #
  # - MacKerell, A.D., Jr., Feig, M., Brooks, C.L., III, Improved Treatment of   #
  #   Protein Backbone Conformational in Empirical Force Fields, J. Am. Chem.    #
  #   Soc. 126(2004): 698-699.                                                   #
  # - MacKerell, A.D., Jr., Feig, M., Brooks, C.L., III, Extending the Treatment #
  #   of Backbone Energetics in Protein Force Fields: Limitations of Gas-Phase   #
  #   Quantum Mechnacis in Reproducing Protein Conformational Distributions in   #
  #   Molecular Dynamics Simulations, J. Comput. Chem. 25(2004): 1400-1415.      #
  # ---------------------------------------------------------------------------- #

  sub CharmmCmap
  {
	print "\nINITIATING CHARMM CMAP SUBROUTINE...\n\n";

	# Reread and analyse $project.data
	my @raw_data;
	open(LAMMPS, "< $project.data") or
	die "\"sub CharmmCmap()\" cannot open \"$project.data!\n";
	print "Analyzing \"$project.data\"...\n\n";
	@raw_data = <LAMMPS>;
	close(LAMMPS);

	# Locate and extract the sections "Masses" and "Atoms"
	my $line_number = 0;
	# Header infos, 0 by default
	my $natom_types = 0;
	my $natom_number = 0;
	my $ndihedral_number = 0;
	my $temp_string;
	# splice points, 0 by default
	my $splice_onset_masses = 0;
	my $splice_onset_atoms = 0;
	my $splice_onset_dihedrals = 0;

	foreach my $line (@raw_data) {
		$line_number++;
		chomp($line);
	# Extract useful informations from the header
		if ($line =~ m/atom types/) {
			($natom_types,$temp_string) = split(" ",$line);
			if ($natom_types == 0) {
				die "\nError: Number of atom types is 0!\n";
			}
			print "Total atom types: $natom_types\n";
		}
		if ($line =~ m/atoms/) {
			($natom_number,$temp_string) = split(" ",$line);
			if ($natom_number == 0) {
				die "\nError: Number of atoms is 0!\n";
			}
			print "Total atoms: $natom_number\n";
		}
		if ($line =~ m/dihedrals/) {
			($ndihedral_number,$temp_string) = split(" ",$line);
			if ($ndihedral_number == 0) {
				die "\nError: Number of dihedrals is 0\n";
			}
			print "Total dihedrals: $ndihedral_number\n";
		}
	# Locate and data from sections "Masses", "Atoms" and "Dihedrals"
		if ($line =~ m/Masses/) {
			$splice_onset_masses = $line_number + 1;
			if ($splice_onset_masses-1 == 0) {
				die "\nError: Can not find the section \"Masses\"\n";
			}
			print "Section \"Masses\" found: line $splice_onset_masses\n";
		}
		if ($line =~ m/Atoms/) {
			$splice_onset_atoms = $line_number +1;
			if ($splice_onset_atoms-1 == 0) {
				die "\nError: Can not find the section \"Atoms\"\n";
			}
			print "Section \"Atoms\" found: line $splice_onset_atoms\n";
		}
		if ($line =~ m/Dihedrals/) {
			$splice_onset_dihedrals = $line_number + 1;
			if ($splice_onset_dihedrals-1 == 0) {
				die "\nError: Can not find the section \"Dihedrals\"\n";
			}
			print "Section \"Dihedrals\" found: line $splice_onset_dihedrals\n";
		}
	}

	print "\nGenerating PHI/PSI dihedral pair list...\n\n";

	my @temp1 = @raw_data;
	my @temp2 = @raw_data;
	my @temp3 = @raw_data;

	# Extract the section "Masses", "Atoms" and "Dihedrals"
	my @temp_masses_data = splice(@temp1,$splice_onset_masses,$natom_types);
	my @temp_atoms_data = splice(@temp2,$splice_onset_atoms,$natom_number);
	my @temp_dihedrals_data = splice(@temp3,$splice_onset_dihedrals,$ndihedral_number);
	
	# Store @temp_masses_dat into a matrix
	my @masses_matrix;
	my $atom_type;
	my $mass;
	for (@temp_masses_data) {
		($atom_type, $mass) = split(" ");
		push(@masses_matrix,[$atom_type,$mass]);
	}

	# Store @temp_atoms_data into a matrix
	my @atoms_matrix;
	my $atom_ID;
	my $molecule_ID;
	my $atype;
	my $charge;
	my $atom_x_coor;
	my $atom_y_coor;
	my $atom_z_coor;
	for (@temp_atoms_data) {
		($atom_ID,$molecule_ID,$atype,$charge,$atom_x_coor,$atom_y_coor,$atom_z_coor) = split(" ");
	 	push(@atoms_matrix,
			[$atom_ID,$molecule_ID,$atype,$charge,$atom_x_coor,$atom_y_coor,$atom_z_coor]);
	}

	# Store @temp_dihedrals_data into a matrix
	my @dihedrals_matrix;
	my $dihedral_ID;
	my $dihedtal_type;
	my $dihe_atom1;
	my $dihe_atom2;
	my $dihe_atom3;
	my $dihe_atom4;
	for (@temp_dihedrals_data) {
		($dihedral_ID,$dihedral_type,$dihe_atom1,$dihe_atom2,$dihe_atom3,$dihe_atom4) = split(" ");
		push(@dihedrals_matrix,
			[$dihedral_ID,$dihedral_type,$dihe_atom1,$dihe_atom2,$dihe_atom3,$dihe_atom4]);
	}
	
	# Find out and extract the peptide backbone dihedrals
	#
	# Definitions of peptide backbone dihedrals
	#
	# For dihedral angle PHI: C--N--CA--C
	# For dihedral angle PSI: N--CA--C--N
	#
	# ---------------------------------------------------------
	#  atom  |  mass  |partial charge|     amino-acid
	# ---------------------------------------------------------
	#   C    | 12.011 |     0.51     | all except GLY and PRO
	#   N    | 14.007 |    -0.29     | PRO
	#   N    | 14.007 |    -0.47     | all except PRO
	#   CA   | 12.011 |     0.07     | all except GLY and PRO
	#   CA   | 12.011 |    -0.02     | GLY
	#   CA   | 12.011 |     0.02     | PRO
	# ---------------------------------------------------------
	#
	#  Peptide backbone
	#        ...
	#         /
	#      O=C
	#         \
	#          N-H
	#         / -----> PHI (C-N-CA-C)
	#      H-CA-R
	#         \ -----> PSI (N-CA-C-N)
	#          C=O
	#         /
	#      H-N
	#         \
	#         ...
	#
	# Criteria to be a PHI/PSI dihedral pair:
	# 1. Atoms have to match with the mass/charge constellations as
	#    defined above.
	# 2. The atoms N--CA--C needs to be covalently bonded with each
	#    other.

	# Find which types do C, N and CA correspond to and store them
	# in lists
	my $mass_carbon = 12.011;
	my $mass_nitrogen = 14.007;

	my @carbon_list;
	my @nitrogen_list;
	my $carbon_counter = 0;
	my $nitrogen_counter = 0;

	for (my $i = 0; $i < $natom_types; $i++) {
		if (${$masses_matrix[$i]}[1] == $mass_carbon) {
			push(@carbon_list,${$masses_matrix[$i]}[0]);
			$carbon_counter++;
		}
		if (${$masses_matrix[$i]}[1] == $mass_nitrogen) {
			push(@nitrogen_list,${$masses_matrix[$i]}[0]);
			$nitrogen_counter++;
		}
	}
	# Quit if no carbons or nitrogens
	if ($carbon_counter == 0 or $nitrogen_counter == 0) {
		if ($carbon_counter == 0) {
			print "No carbon atoms exist in the system\n";
		}
		if ($nitrogen_counter == 0) {
			print "No nitrogen atoms exist in the system\n";
		}
		print "CMAP usage impossible\n";
		return;
	}

	print "Carbon atom type/s: @carbon_list\n";
	print "Nitrogen atom type/s: @nitrogen_list\n";

	# Determine the atom types of C, CA and N

	# Charges of the backbone atoms
	my $charge_C = 0.51;
	my $charge_CA = 0.07;
	my $charge_N = -0.47;

	# Special setting for PRO
	my $charge_N_PRO = -0.29;
	my $charge_CA_PRO = 0.02;

	# Special setting for GLY
	my $charge_CA_GLY = -0.02;

	# Peptide backbone atom types
	my $C_type;
	my $CA_type;
	my $CA_GLY_type;
	my $CA_PRO_type;
	my $N_type;
	my $N_PRO_type;

	my $C_counter = 0;
	my $CA_counter = 0;
	my $CA_GLY_counter = 0;
	my $CA_PRO_counter = 0;
	my $N_counter = 0;
	my $N_PRO_counter = 0;

	my $C_flag = 0;

	for (my $i = 0; $i <= $natom_number; $i++) {
		my $cur_type = ${$atoms_matrix[$i]}[2];
		my $cur_charge = ${$atoms_matrix[$i]}[3];
		for (my $j = 0; $j <= $#carbon_list; $j++) {
			if ($cur_type == $carbon_list[$j]) {
				$C_flag = 1;
				if ($cur_charge == $charge_C) {
					$C_type = $cur_type;
					$C_counter++;
				}
				if ($cur_charge == $charge_CA) {
					$CA_type = $cur_type;
					$CA_counter++;
				}
				if ($cur_charge == $charge_CA_GLY) {
					$CA_GLY_type = $cur_type;
					$CA_GLY_counter++;
				}
				if ($cur_charge == $charge_CA_PRO) {
					$CA_PRO_type = $cur_type;
					$CA_PRO_counter++;
				}
			}
		}
		if ($C_flag == 0) {
			for (my $k = 0; $k <= $#nitrogen_list; $k++) {
				if ($cur_type == $nitrogen_list[$k]) {
					if ($cur_charge == $charge_N) {
						$N_type = $cur_type;
						$N_counter++;
					}
					if ($cur_charge == $charge_N_PRO) {
						$N_PRO_type = $cur_type;
						$N_PRO_counter++;
					}
				}
			}
		}
		$C_flag = 0;
	}

	# Quit if one of the atom types dosen't exist
	if ( $C_counter == 0 or
	    ($CA_counter == 0 and $CA_GLY_counter == 0 and $CA_PRO_counter == 0) or
	    ($N_counter == 0 and $N_PRO_counter == 0) ) {
		if ($C_counter == 0) {
			print "\nCannot find the peptide backbone C atom type\n";
		}
		if ($CA_counter == 0 and $CA_GLY_counter == 0 and $CA_PRO_counter == 0) {
			print "\nCannot find the peptide backbone C-alpha atom type\n";
		}
		if ($N_counter == 0 and $N_PRO_counter == 0) {
			print "\nCannot find the peptide backbone N atom type\n";
		}
		print "CMAP usage impossible\n";
		return;
	}

	print "Peptide backbone carbon type: $C_type\n";
	print "Alpha-carbon type: $CA_type\n" if ($CA_counter > 0);
	print "Alpha-carbon type (GLY): $CA_GLY_type\n" if ($CA_GLY_counter > 0);
	print "Alpha-carbon type (PRO): $CA_PRO_type\n" if ($CA_PRO_counter > 0);
	print "Peptide backbone nitrogen type: $N_type\n" if ($N_counter >0);
	print "Peptide backbone nitrogen type (PRO): $N_PRO_type\n" if ($N_PRO_counter > 0);

	# Loop through the dihedral list to find the PHI- and PSI-dihedrals
	my @PHI_dihedrals;
	my @PSI_dihedrals;
	my $PHI_counter = 0;
	my $PSI_counter = 0;

	for (my $i = 0; $i < $ndihedral_number; $i++) {
		my $cur_dihe_ID = ${dihedrals_matrix[$i]}[0];
		my $cur_atom1_type = ${atoms_matrix[${dihedrals_matrix[$i]}[2]-1]}[2];
		my $cur_atom2_type = ${atoms_matrix[${dihedrals_matrix[$i]}[3]-1]}[2];
		my $cur_atom3_type = ${atoms_matrix[${dihedrals_matrix[$i]}[4]-1]}[2];
		my $cur_atom4_type = ${atoms_matrix[${dihedrals_matrix[$i]}[5]-1]}[2];

		next if (${dihedrals_matrix[$i]}[2] == ${dihedrals_matrix[$i-1]}[2] and
		         ${dihedrals_matrix[$i]}[3] == ${dihedrals_matrix[$i-1]}[3] and
	                 ${dihedrals_matrix[$i]}[4] == ${dihedrals_matrix[$i-1]}[4] and
                         ${dihedrals_matrix[$i]}[5] == ${dihedrals_matrix[$i-1]}[5]);

		# Determine PHI-dihedrals; If C-CA-N-C or C-N-CA-C, then save it in a list
		if ($cur_atom1_type == $C_type and $cur_atom4_type == $C_type) {
			if ( ( ($cur_atom2_type == $CA_type or
			        $cur_atom2_type == $CA_GLY_type or
				$cur_atom2_type == $CA_PRO_type) and
			       ($cur_atom3_type == $N_type or
			        $cur_atom3_type == $N_PRO_type) ) or
			     ( ($cur_atom3_type == $CA_type or
				$cur_atom3_type == $CA_GLY_type or
			        $cur_atom3_type == $CA_PRO_type) and
			       ($cur_atom2_type == $N_type or
				$cur_atom2_type == $N_PRO_type) ) ) {
				push (@PHI_dihedrals,$cur_dihe_ID);
				$PHI_counter++;
			}
		}

		# Determin PSI-dihedrals; If N-CA-C-N or N-C-CA-N (N can be both normal N or N proline),
		# then save it in a list
		if ( ($cur_atom1_type == $N_type and $cur_atom4_type == $N_type) or
		     ($cur_atom4_type == $N_PRO_type and $cur_atom1_type == $N_PRO_type) or
	     	     ($cur_atom1_type == $N_type and $cur_atom4_type == $N_PRO_type) or
                     ($cur_atom4_type == $N_type and $cur_atom1_type == $N_PRO_type) ) {
			if ( ( ($cur_atom2_type == $CA_type or
			        $cur_atom2_type == $CA_GLY_type or
		                $cur_atom2_type == $CA_PRO_type) and
			        $cur_atom3_type == $C_type) or
			     ( ($cur_atom3_type == $CA_type or
				$cur_atom3_type == $CA_GLY_type or
			        $cur_atom3_type == $CA_PRO_type) and
				$cur_atom2_type == $C_type) ) {
				push (@PSI_dihedrals,$cur_dihe_ID);
				$PSI_counter++;
			}
		}
	}

	# Quit if no PHI or PSI dihedrals
	if ($PHI_counter == 0 or $PSI_counter ==0) {
		if ($PHI_counter == 0) {
			print "Can not find the PHI backbone dihedrals\n";
		}
		if ($PSI_counter ==0) {
			print "Can not find the PSI backbone dihedrals\n";
		}
		print "CMAP usage impossible\n";
		return;
	}

	# Construct the PHI/PSI diheral pair list
	#
	# The algorithm:
	#         _____
	#        |     |
	#     1--2--3--4      PHI-dihedral
	#     4--3--2--1
	#   --C--N-CA--C--N-- Peptide backbone
	#        1--2--3--4
	#        4--3--2--1   PSI-dihedral
	#        |_____|
	#
	# For a certain PHI dihedral, following conditions have to be met:
	#
	#        PHI         PSI
	#  If (2--3--4) = (1--2--3)
	#  or
	#  if (2--3--4) = (4--3--2)
	#  or
	#  if (3--2--1) = (1--2--3)
	#  or
	#  if (3--2--1) = (4--3--2),
	#
	# then these 2 dihedrals are a PHI/PSI pair. If a pair is found, the
	# dihedral IDs will be stored in "@PHI_PSI_matrix".
	
	my @PHI_PSI_matrix;
	my $crossterm_CA_charge;
	my $crossterm_type;
	my $crossterm_counter = 0;
	my $crossterm_type1_flag = 0;
        my $crossterm_type2_flag = 0;
	my $crossterm_type3_flag = 0;
	my $crossterm_type4_flag = 0;
	my $crossterm_type5_flag = 0;
	my $crossterm_type6_flag = 0;

	for (my $i = 0; $i <= $#PHI_dihedrals; $i++) {
		my $cur_PHI_dihe = ${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[0];
		my $phi1 = ${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[2];
		my $phi2 = ${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[3];
		my $phi3 = ${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[4];
		my $phi4 = ${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[5];
		for (my $j = 0; $j <= $#PSI_dihedrals; $j++) {
			my $cur_PSI_dihe = ${dihedrals_matrix[$PSI_dihedrals[$j]-1]}[0];
			my $psi1 = ${dihedrals_matrix[$PSI_dihedrals[$j]-1]}[2];
			my $psi2 = ${dihedrals_matrix[$PSI_dihedrals[$j]-1]}[3];
			my $psi3 = ${dihedrals_matrix[$PSI_dihedrals[$j]-1]}[4];
			my $psi4 = ${dihedrals_matrix[$PSI_dihedrals[$j]-1]}[5];
			if ( ($phi2 == $psi1 and $phi3 == $psi2 and $phi4 == $psi3) or
			     ($phi2 == $psi4 and $phi3 == $psi3 and $phi4 == $psi2) or
			     ($phi3 == $psi1 and $phi2 == $psi2 and $phi1 == $psi3) or
			     ($phi3 == $psi4 and $phi2 == $psi3 and $phi1 == $psi2) ) {

			     # Find out to which amino acid the cross-term belongs

			     if ($phi3 == $psi2 or $phi3 == $psi3) {
				     $crossterm_CA_charge = ${atoms_matrix[${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[4]-1]}[3];
			     }
			     if ($phi2 == $psi2 or $phi2 == $psi3) {
				     $crossterm_CA_charge = ${atoms_matrix[${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[3]-1]}[3];
			     }

			     # Def. the type of the crossterm per cmap.data file; If C_alpha of the crossterm is
			     # - ALA type, then $crossterm_type = 1;
                             # - ALA-PRO (ALA is the current AA), then $crossterm_type = 2;
			     # - PRO type, then $crossterm_type = 3;
                             # - PRO-PRO (First PRO is the current AA), then $crossterm_type = 4;
			     # - GLY type, then $crossterm_type = 5;
                             # - GLY-PRO (GLY is the current AA), then $crossterm_type = 6;

			     if ($crossterm_CA_charge == $charge_CA) { $crossterm_type = 1; $crossterm_type1_flag = 1; }
			     if ($crossterm_CA_charge == $charge_CA_GLY) { $crossterm_type = 5; $crossterm_type5_flag = 1; }
			     if ($crossterm_CA_charge == $charge_CA_PRO) {
					$crossterm_type = 3; $crossterm_type3_flag = 1;
					# Checking the last crossterm, re-assign the last crossterm type if needed
					if ($crossterm_counter-1 >= 0 and $PHI_PSI_matrix[$crossterm_counter-1][0] == 1) {
						$PHI_PSI_matrix[$crossterm_counter-1][0] = 2;
						$crossterm_type2_flag = 1;
					}
					if ($crossterm_counter-1 >= 0 and $PHI_PSI_matrix[$crossterm_counter-1][0] == 3) {
                                                $PHI_PSI_matrix[$crossterm_counter-1][0] = 4;
                                                $crossterm_type4_flag = 1;
                                        }
					if ($crossterm_counter-1 >= 0 and $PHI_PSI_matrix[$crossterm_counter-1][0] == 5) {
                                                $PHI_PSI_matrix[$crossterm_counter-1][0] = 6;
                                                $crossterm_type6_flag = 1;
                                        }
			     }	
			     push(@PHI_PSI_matrix,[$crossterm_type,$phi1,$phi2,$phi3,$phi4,$psi4]);
			     $crossterm_counter++;

			     $crossterm_CA_charge = 0;
			     $crossterm_type = 0;
		     }
		}
	}

	# Check whether the amino acid at the C-terminus is a PRO or not. If yes, the type of the last crossterm
	# should be set to its X-PRO form  instead of X, where X is ALA, PRO, or GLY.  X-PRO form = X type + 1.

	my @pdb_data;
	open(PDB,"< $project.pdb")
		or die "WARNING: Cannot open file \"$project.pdb\"! (required if the -cmap option is used)\n";
	@pdb_data = <PDB>;
	close(PDB);

        my @ter_line;
        my $ter_AA_type = 0;
        my $ter_flag = 0;
	foreach $line (@pdb_data) {
		if ($line =~ m/TER/) {
		      @ter_line = split(" ",$line);
                      $ter_AA_type = $ter_line[2]; 	
                      print "Terminal amino acid type is: $ter_AA_type\n";
                      $ter_flag = 1;
	        }
        }
        if ($ter_flag == 0) {
                print "\n*** ERROR IN THE PDB FILE: ***\n";
                print "In order for the CMAP section to be generated, the pdb file must \n";
                print "identify the C-terminus amino acid in the file with 'TER'. \n";
                print "This line is missing from the pdb file that was used.\n";
                print "To correct this problem, open the pdb file in an editor,\n";
                print "find the last atom of the last amino acid residue in the peptide\n";
                print "chain and insert the following line immediately after that atom:\n";
                print "                  'TER   <#1>   <RES>   <#2>' \n";
                print "where '<#1> is the next atom number, <RES> is the three letter amino\n";
                print "acid abbreviation for that amino acid, and <#2> is the molecule number\n";
                print "of the terminal amino acid residue.\n\n";
                print "For example, if the last atom of the last amino acid in the peptide\n";
                print "sequence is listed in the pdb file as:\n\n";
                print "  'ATOM   853  O  GLU  P  56  12.089 -1.695 -6.543 1.00 1.03  PROA'\n\n";
                print "you would insert the following line after it:\n\n";
                print "  'TER    854     GLU     56'\n\n";
                print "If any additional atoms are listed in the pdb file (e.g., water, ions)\n";
                print "after this terminal amino acid residue, their atom numbers and\n";
                print "molecule numbers must be incremented by 1 to account for the new line\n";
                print "that was inserted.\n\n";
                die "Error: No terminating atom designated in pdb file!  See above note to correct problem.\n\n";
        }

        if ($ter_AA_type eq PRO) {
                $PHI_PSI_matrix[$crossterm_counter-1][0] = $PHI_PSI_matrix[$crossterm_counter-1][0]+1;
        }

	# Print out PHI/PSI diheral pair list
	my $pair_counter = 0;
        # Don't presently use $ncrosstermtypes but have this available if wish to print it out
	my $ncrosstermtypes = $crossterm_type1_flag + $crossterm_type2_flag + $crossterm_type3_flag +
                              $crossterm_type4_flag + $crossterm_type5_flag + $crossterm_type6_flag;
	print "\nWriting \"$project.data\" with section \"CMAP crossterms\" added at the end.\n";
	
	# Writing the new lammps data file
	open(REWRITE,"> $project.data")
		or die "Cannot write file \"$project.data\"!\n";
	foreach $line (@raw_data) {
		printf(REWRITE "$line\n");
		if ($line =~ m/impropers/) {
			printf(REWRITE "%12d  crossterms\n", $crossterm_counter);
		}
	}
	printf(REWRITE "CMAP\n\n");

        my $ref_line;
        my $column;
	foreach $ref_line (@PHI_PSI_matrix) {
		$pair_counter++;
		printf(REWRITE "%8d",$pair_counter);
		foreach $column (@$ref_line) {
			printf(REWRITE " %7d",$column);
		}
		printf(REWRITE "\n");
	}
	close(REWRITE);
	print "\nDone!\n\n";
}

# main

  Initialize();
  WriteData();
  WriteLAMMPSInput();
  printf("Info: conversion complete\n\n") if ($info);
  CharmmCmap() if ($cmap);
