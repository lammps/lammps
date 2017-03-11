#!/usr/bin/perl
#
#  program:	lammps2pdb.pl
#  author:	Pieter J. in 't Veld
#  		pjintve@sandia.gov, veld@verizon.net
#  date:	September 3-4, 2004, April 23-25, 2005.
#  purpose:	Translation of lammps output into pdb format for use in
#  		conjunction with vmd
#
#  Notes:	Copyright by author for Sandia National Laboratories
#    20040903	Conception date of v1.0: rudimentary script for collagen
#    		project.
#    20050423	Conception date of v2.0:
#    		- changed data access through indexing data directly on disk;
#    		- added all command line options
#    20050425	Corrected focussing to use a target molecule's moment of
#    		inertia to determine its principal axis: for computational
#    		ease, a 3D orientational vector is calculated from two 2D
#    		projections
#    20050701	Changed create_names() to loosen mass recognition criterion and
#    		corrected created_atom_ids() to properly deal with the end of
#    		the data stream.

# subroutines

  sub test
  {
    my $name		= shift(@_);

    printf("Error: file %s not found\n", $name) if (!scalar(stat($name)));
    return !scalar(stat($name));
  }


  sub initialize
  {
    my $k		= 0;
    my @options		= ("-help", "-nstart", "-dn", "-cut", "-repair",
			   "-units", "-quiet", "-nodetect", "-data", "-pbc",
			   "-focus", "-center", "-exclude");
    my @remarks		= ("display this message",
			   "starting frame [-1]",
			   "frames to skip when creating multiple frames [0]",
			   "cut bonds crossing over box edge [off]",
			   "repair bonds [off]",
			   "dump file entries have units [off]",
			   "turn display information off",
			   "turn auto-detection of masses off",
			   "use data file other than project.data",
			   "apply periodic boundary to molecules []",
			   "molecules to focus on []",
			   "molecules to use as center of mass []",
			   "exclude molecules from output []");
    my @notes		= (
      "* Multiple frames are processed when dn > 0.",
      "* Only the last frame is converted when nstart < 0.",
      "* Focussing superceedes centering.",
      "* Focussing uses the moment of inertia to determine the molecules'",
      "  principal axis; the molecules are rotated and translated to the lab",
      "  frame, using the focus molecules as reference.",
      "* Expected files in current directory: project.data, project.dump",
      "* Generated output files: project_trj.psf, project_trj.pdb",
      "* Uses project_ctrl.psf if available");


    $program		= "lammps2pdb";
    $version		= "2.2.5";
    $year		= "2007";
    $cut		= 0;
    $repair		= 0;
    $units		= 0;
    $info		= 1;
    $nstart		= -1;
    $dn			= 0;
    $detect		= 1;

    foreach (@ARGV)
    {
      if (substr($_, 0, 1) eq "-")
      {
	my $k		= 0;
	my @arg		= split("=");
	my $switch	= (($arg[1] eq "")||($arg[1] eq "on")||($arg[1]!=0));
	$arg[0]		= lc($arg[0]);
	foreach (@options)
	{
	  last if ($arg[0] eq substr($_, 0, length($arg[0])));
	  ++$k;
	}
	$help		= 1 if (!$k--);
	$nstart		= $arg[1] if (!$k--);
	$dn		= $arg[1] if (!$k--);
	$cut		= $switch if (!$k--);
	$repair		= $switch if (!$k--);
	$units		= $switch if (!$k--);
	$info		= $switch ? 0 : 1 if (!$k--);
	$detect		= $switch ? 0 : 1 if (!$k--);
	$data_name	= $arg[1] if (!$k--);
	if (!$k--) {
	  if ($switch) { $pbc{ALL} = 1; }
	  else { foreach (split(",",$arg[1])) { $pbc{uc($_)} = 1; }}}
	if (!$k--) { foreach (split(",",$arg[1])) { $focus{uc($_)} = uc($_);}}
	if (!$k--) { foreach (split(",",$arg[1])) { $center{uc($_)} = uc($_);}}
	if (!$k--) { foreach (split(",",$arg[1])) { $exclude{uc($_)}=uc($_);}}
      }
      else { $project	= $_ if (!$k++); }
    }
    if (($k<1)||$help)
    {
      printf("%s v%s (c)%s by Pieter J. in \'t Veld for SNL\n\n",
        $program, $version, $year);
      printf("Usage:\n  %s.pl [-option[=#] ..] project\n\n", $program);
      printf("Options:\n");
      for (my $i=0; $i<scalar(@options); ++$i)
      {
	printf("  %-10.10s %s\n", $options[$i], $remarks[$i]);
      }
      printf("\nNotes:\n") if (scalar(@notes));
      foreach(@notes) { printf("  %s\n", $_); };
      printf("\n");
      exit(-1);
    }
    printf("%s v%s (c)%s\n\n", $program, $version, $year) if ($info);

    $data_name		= $project.".data" if ($data_name eq "");
    $traject_name	= $project.".dump";
    $pdb_name		= $project."_trj.pdb";
    $psf_name		= $project."_trj.psf";
    $psf_ctrl_name	= $project."_ctrl.psf";
    $psf_ctrl		= scalar(stat($psf_ctrl_name));

    $traject_flag	= scalar(stat($traject_name));
    my $flag		= test($data_name);
    printf("Conversion aborted\n\n") if ($flag);
    exit(-1) if ($flag);

    # data input

    create_atom_ids();
    create_bonds();
    open(TRAJECT, "<".$traject_name) if ($traject_flag);
    open(PSF_CTRL, "<".$psf_ctrl_name) if ($psf_ctrl);

    # open output

    open(PSF, ">".$psf_name) if (!$psf_ctrl);
    open(PDB, ">".$pdb_name);

    # align center with focus

    %center		= %focus if (scalar(%focus));
  }


# Vector routines

  sub V_String					# V_String(@a)
  {
    return "{".join(", ", @_)."}";
  }


  sub V_Round
  {
    foreach(@_) { $_ = ($_<0 ? -int(-$_*1e11+0.5) : int($_*1e11+0.5))/1e11; }
    return @_;
  }


  sub V_Add					# V_Add(@a, @b) = @a + @b;
  {
    return (@_[0]+@_[3], @_[1]+@_[4], @_[2]+@_[5]);
  }


  sub V_Subtr					# V_Subtr(@a, @b) = @a - @b;
  {
    return (@_[0]-@_[3], @_[1]-@_[4], @_[2]-@_[5]);
  }


  sub V_Dot					# V_Dot(@a, @b) = @a . @b;
  {
    return (@_[0]*@_[3]+@_[1]*@_[4]+@_[2]*@_[5]);
  }


  sub V_Cross					# V_Cross(@a, @b) = @a x @b;
  {
    return (@_[1]*@_[5]-@_[2]*@_[4], @_[2]*@_[3]-@_[0]*@_[5],
      @_[0]*@_[4]-@_[1]*@_[3]);
  }


  sub V_Mult					# V_Mult($f, @a) = $f * @a;
  {
    return (@_[0]*@_[1], @_[0]*@_[2], @_[0]*@_[3]);
  }


  sub V_Norm					# V_Norm(@a) = @a/|@a|;
  {
    return V_Mult(1/sqrt(V_Dot(@_[0,1,2],@_[0,1,2])), @_[0,1,2]);
  }


  sub pbc					# periodic -0.5*L <= x < 0.5*L
  {
    my $x		= @_[0]/@_[1]+0.5;

    return @_[0]-@_[1]*($x<0 ? int($x)-1 : int($x));
  }


  sub V_PBC					# V_PBC(@a, @l)
  {
    return (pbc(@_[0], @_[3]), pbc(@_[1], @_[4]), pbc(@_[2], @_[5]));
  }


  sub V_Cmp					# V_Cmp(abs(@a), abs(@b))
  {
    return -1 if (abs($_[0])<abs($_[3]));
    return  1 if (abs($_[0])>abs($_[3]));
    return -1 if (abs($_[1])<abs($_[4]));
    return  1 if (abs($_[1])>abs($_[4]));
    return -1 if (abs($_[2])<abs($_[5]));
    return  1 if (abs($_[2])>abs($_[5]));
    return 0;
  }


  sub V_Sort					# sort on descending absolute
  {						# value
    my @v		= @_;

    for (my $i=0; $i<scalar(@v)-3; $i+=3)
    {
      for (my $j=$i+3; $j<scalar(@v); $j+=3)
      {
	my @u		= @v[$i, $i+1, $i+2];
	next if (V_Cmp(@u, @v[$j, $j+1, $j+2])>=0);
	@v[$i,$i+1,$i+2]= @v[$j,$j+1,$j+2];
	@v[$j,$j+1,$j+2]= @u;
      }
    }
    return @v;
  }


# Matrix routines

  sub M_String					# M_String(@A)
  {
    my @M;

    for (my $i=0; $i<3; ++$i) { push(@M, V_String(splice(@_, 0, 3))); }
    return "{".join(", ", @M)."}";
  }


  sub M_Round
  {
    return V_Round(@_);
  }


  sub M_Transpose				# M_Transpose(@A) = (@A)^t;
  {
    return (@_[0], @_[3], @_[6], @_[1], @_[4], @_[7], @_[2], @_[5], @_[8]);
  }


  sub M_Dot					# M_Dot(@A, @B) = @A . @B;
  {
    return (
      V_Dot(@_[0,1,2], @_[ 9,12,15]), V_Dot(@_[0,1,2], @_[10,13,16]),
	V_Dot(@_[0,1,2], @_[11,14,17]),
      V_Dot(@_[3,4,5], @_[ 9,12,15]), V_Dot(@_[3,4,5], @_[10,13,16]),
	V_Dot(@_[3,4,5], @_[11,14,17]),
      V_Dot(@_[6,7,8], @_[ 9,12,15]), V_Dot(@_[6,7,8], @_[10,13,16]),
	V_Dot(@_[6,7,8], @_[11,14,17]));
  }


  sub M_Det					# M_Det(@A) = | @A |
  {
    return V_Dot(@_[0,1,2], V_Cross(@_[3,4,5], @_[6,7,8]));
  }


  sub M_Mult					# M_Mult($a, @A) = $a * @A
  {
    return (
      @_[0]*@_[1], @_[0]*@_[2], @_[0]*@_[3],
      @_[0]*@_[4], @_[0]*@_[5], @_[0]*@_[6],
      @_[0]*@_[7], @_[0]*@_[8], @_[0]*@_[9]);
  }


  sub M_Unit { return (1,0,0, 0,1,0, 0,0,1); }

  sub PI { return 4*atan2(1,1); }

  sub M_Rotate					# M_Rotate($n, $alpha) = @R_$n;
  {						# vmd convention
    my $n		= shift(@_);
    my $alpha		= shift(@_)*PI()/180;
    my $cos		= cos($alpha);
    my $sin		= sin($alpha);

    $cos		= 0 if (abs($cos)<1e-16);
    $sin		= 0 if (abs($sin)<1e-16);
    return (1,0,0, 0,$cos,-$sin, 0,$sin,$cos) if ($n==0); # around x-axis
    return ($cos,0,$sin, 0,1,0, -$sin,0,$cos) if ($n==1); # around y-axis
    return ($cos,-$sin,0, $sin,$cos,0, 0,0,1) if ($n==2); # around z-axis
    return M_Unit();
  }


  sub M_RotateNormal				# returns R.(1,0,0) = @a/|@a|;
  {
    my @R		= M_Unit();
    my @n		= V_Mult(1.0/sqrt(V_Dot(@_[0,1,2], @_)), @_);
    my $sina		= $n[1]<0 ? -sqrt($n[1]*$n[1]+$n[2]*$n[2]) :
				     sqrt($n[1]*$n[1]+$n[2]*$n[2]);

    if ($sina)
    {
      my $cosa		= $n[0];
      my $cosb		= $n[1]/$sina;
      my $sinb		= $n[2]/$sina;
      @R		= (
	$cosa      ,-$sina      ,      0,
	$cosb*$sina, $cosb*$cosa, -$sinb,
	$sinb*$sina, $sinb*$cosa,  $cosb);
    }
    return @R;
  }


  sub MV_Dot					# MV_Dot(@A, @b) = @A . @b;
  {
    return (V_Dot(@_[0,1,2], @_[9,10,11]), V_Dot(@_[3,4,5], @_[9,10,11]),
      V_Dot(@_[6,7,8], @_[9,10,11]));
  }


  sub M_Sweep
  {
    my @M		= @_;

    for (my $i=0; $i<3; ++$i)
    {
      @M		= V_Sort(@M);
      my @v		= splice(@M, 3*$i, 3);
      if (abs($v[$i])>1e-10)
      {
	@v		= V_Mult(1/$v[$i], @v);
	@M[0,1,2]	= V_Subtr(@M[0,1,2], V_Mult($M[$i], @v));
	@M[3,4,5]	= V_Subtr(@M[3,4,5], V_Mult($M[$i+3], @v));
      }
      @M		= (@v, @M) if ($i==0);
      @M		= (@M[0,1,2], @v, @M[3,4,5]) if ($i==1);
      @M		= (@M, @v) if ($i==2);
    }
    return V_Round(@M);
  }


# Complex routines

  sub C_Add					# z = z1 + z2
  {
    return (@_[0]+@_[2], @_[1]+@_[3]);
  }


  sub C_Subtr					# z = z1 - z2
  {
    return (@_[0]-@_[2], @_[1]-@_[3]);
  }


  sub C_Mult					# z = z1 * z2
  {
    return (@_[0]*@_[2]-@_[1]*@_[3], @_[0]*@_[3]+@_[2]*@_[1]);
  }


  sub C_Div					# z = z1 / z2
  {
    my $num		= @_[2]*@_[2]+@_[3]*@_[3];

    return ((@_[0]*@_[2]+@_[1]*@_[3])/$num, (@_[1]*@_[2]-@_[0]*@_[3])/$num);
  }


  sub C_Pow					# z = z1 ^ a
  {
    my $r		= (sqrt(@_[0]*@_[0]+@_[1]*@_[1]))**@_[2];
    my $a		= @_[2]*atan2(@_[1], @_[0]);

    return ($r*cos($a), $r*sin($a));
  }


  sub C_Correct
  {
    return (abs(@_[0])<1e-14 ? 0 : @_[0], abs(@_[1])<1e-14 ? 0 : @_[1]);
  }


  sub C_Zero
  {
    return (abs(@_[0])<1e-8 ? 0 : @_[0], abs(@_[1])<1e-8 ? 0 : @_[1]);
  }


  sub C_Conj
  {
    return (@_[0], -@_[1]);
  }


  sub C_String
  {
    return @_[0]." + ".@_[1]."i";
  }


# analytical roots

  sub R_Sort
  {
    my $n		= scalar(@_);

    for (my $i=0; $i<$n-2; $i+=2)
    {
      for (my $j=$i+2; $j<$n; $j+=2)
      {
	if (@_[$j]<@_[$i]) {
	  my @t = @_[$i,$i+1]; @_[$i,$i+1] = @_[$j,$j+1]; @_[$j,$j+1] = @t; }
	else { if ((@_[$j]==@_[$i])&&(@_[$j+1]<@_[$i+1])) {
	  my @t = @_[$i,$i+1]; @_[$i,$i+1] = @_[$j,$j+1]; @_[$j,$j+1] = @t; } }
      }
    }
    return @_;
  }


  sub R_First
  {
    return (0, 0) if (abs(@_[1])<1e-14);
    return (-@_[0]/@_[1], 0);
  }


  sub R_Second
  {
    return R_First(@_) if (abs(@_[2]<1e-14));
    my @z		= (-@_[1]/@_[2]/2.0, 0);
    my @root		= C_Correct(C_Pow(($z[0]*$z[0]-@_[0]/@_[2], 0), 1/2));
    return (C_Zero(C_Subtr(@z, @root)), C_Zero(C_Add(@z, @root)));
  }


  sub R_Third
  {
    return R_Second(@_) if (abs(@_[3])<1e-14);
    my $c3		= 3*@_[3];
    my $f1		= @_[1]*$c3-@_[2]*@_[2];
    my $f2		= 0.5*@_[2]*(3*$f1+@_[2]*@_[2])-1.5*@_[0]*$c3*$c3;
    my $f3		= -@_[2]/$c3;
    my @A		= (0, 0);
    my @B		= (0, 0);
    my @z1		= (0.5, 0.5*sqrt(3));
    my @z2		= C_Conj(@z1);

    if (abs($f1)<1e-3) { 				# limit f1 -> 0
      @A		= ($f2<0 ? abs(2*$f2)**(1/3) : 0, 0); }
    else {
      if (abs($f2)<1e-14) {				# limit f2 -> 0
	my $f		= sqrt(abs($f1))/$c3;
	@A		= $f1<0 ? (-$f*$z1[1], 0.5*$f) : ($f, 0);
	@B		= $f1<0 ? (-$A[0], $A[1]) : ($f, 0); }
      else {
	@B		= C_Pow(C_Add(($f2, 0),
			    C_Pow(($f1*$f1*$f1+$f2*$f2, 0), 1/2)), 1/3);
	@A		= C_Div(($f1/$c3, 0), @B);
	@B		= ($B[0]/$c3, $B[1]/$c3); } }
    return R_Sort(
      C_Zero(C_Add(($f3, 0), C_Subtr(C_Mult(@z1, @A), C_Mult(@z2, @B)))),
      C_Zero(C_Add(($f3, 0), C_Subtr(C_Mult(@z2, @A), C_Mult(@z1, @B)))),
      C_Zero(C_Add(($f3, 0), C_Subtr(@B, @A))));
  }


# general read file operators

  sub advance_read
  {
    my $input		= shift;
    my $dlines		= shift;
    my $read;

    while ($dlines--) { $read = <$input>; }
    return $read;
  }


  sub rewind_read
  {
    my $input		= shift;

    seek($input, 0, SEEK_SET);
  }


# create crossreference tables

  sub create_psf_index					# make an index of id
  {							# locations
    my @psf_ids		= ("!NATOM","!NBOND:","!NTHETA:","!NPHI:","!NIMPHI:");
    my @ids		= (atoms, bonds, angles, dihedrals, impropers);
    my $k		= 0;
    my %hash;

    printf("Info: creating PSF index\n") if ($info);
    foreach (@psf_ids) { $hash{$_} = shift(@ids); };
    seek(PSF_CTRL, 0, SEEK_SET);
    while (<PSF_CTRL>)
    {
      chop();
      my @tmp		= split(" ");
      my $n		= $hash{$tmp[1]};
      $PSFIndex{$n}	= tell(PSF_CTRL)." ".$tmp[0] if ($n ne "");
    }
  }


  sub psf_goto						# goto $ident in <PSF>
  {
    create_psf_index() if (!scalar(%PSFIndex));
    my $id		= shift(@_);
    my @n 		= split(" ", $PSFIndex{$id});

    @PSFBuffer		= ();
    if (!scalar(@n))
    {
      printf("Warning: PSF index for $id not found\n");
      seek(PSF_CTRL, 0, SEEK_END);
      return -1;
    }
    seek(PSF_CTRL, $n[0], SEEK_SET);
    return $n[1];
  }


  sub psf_get
  {
    @PSFBuffer		= split(" ", <PSF_CTRL>) if (!scalar(@PSFBuffer));
    return splice(@PSFBuffer, 0, shift(@_));
  }


  sub create_data_index
  {
    my $init		= 3;
    my @tmp;
    my $id;
    my %hash;
    my %size;

    foreach ((masses,atoms,bonds,angles,dihedrals,impropers)) { $hash{$_}=$_; }
    open(DATA, "<".$data_name);
    for (my $i=0; $i<2; ++$i) { my $tmp = <DATA>; }	# skip first two lines
    while ($init&&!eof(DATA))				# interpret header
    {
      @tmp		= split(" ", <DATA>);
      --$init if (!scalar(@tmp));
      foreach (@tmp) { $_ = lc($_); }
      $tmp[1]		= masses if ($tmp[1]." ".$tmp[2] eq "atom types");
      $L[0]		= $tmp[0]." ".($tmp[1]-$tmp[0])
					if (join(" ",@tmp[2,3]) eq "xlo xhi");
      $L[1]		= $tmp[0]." ".($tmp[1]-$tmp[0])
					if (join(" ",@tmp[2,3]) eq "ylo yhi");
      $L[2]		= $tmp[0]." ".($tmp[1]-$tmp[0])
					if (join(" ",@tmp[2,3]) eq "zlo zhi");
      if (($id = $hash{$tmp[1]}) ne "") { $size{$id}	= $tmp[0]; }
    }
    @l			= ();
    for (my $i=0; $i<3; ++$i)
    {
      @tmp		= split(" ", $L[$i]);
      $l[$i]		= $tmp[0];
      $l[$i+3]		= $tmp[1];
    }
    while (!eof(DATA))					# interpret data
    {
      @tmp		= split(" ", <DATA>);
      my $skip		= <DATA> if (($id = $hash{lc($tmp[0])}) ne "");
      $DATAIndex{$id}	= tell(DATA)." ".$size{$id} if ($id ne "");
    }
  }


  sub goto_data
  {
    create_data_index() if (!scalar(%DATAIndex));
    my $id		= shift(@_);
    my @n		= split(" ", $DATAIndex{$id});

    if (!scalar(@n))
    {
      printf("Warning: DATA index for $id not found\n");
      seek(DATA, 0, SEEK_END);
      return -1;
    }
    seek(DATA, $n[0], SEEK_SET);
    return $n[1];
  }


# create atom and residue identifiers

  sub create_names
  {
    return if (scalar(@names));
    my $n		= goto_data(masses);
    my @mass		= (1, 12, 14, 16, 32.1, 4, 20.2, 40.1, 65.4,
			   55.8, 31, 35.5, 23, 24.3, 39.1, 28.1);
    my @name 		= ("H", "C", "N", "O", "S", "HE", "NE", "CA",
			   "ZN", "FE", "P", "CL", "NA", "MG", "K", "SI");
    my @unknown		= ("X", "Y", "Z");
    my @letter		= ("X", "Y", "Z", "P", "Q", "R", "S", "A", "B", "C");
    my $k		= 0;
    my %atom;

    $names[0]		= "";
    foreach (@mass) { $atom{$_} = shift(@name); }
    for (my $i=1; $i<=$n; ++$i)
    {
      my @tmp		= split(" ", <DATA>);
      my $j		= $tmp[0];
      $tmp[1]		= int($tmp[1]*10+0.5)/10;
      if ((($names[$j] = $atom{$masses[$j] = $tmp[1]}) eq "")||!$detect)
      {
	$names[$j]	= $letter[$k].$unknown[0];
	if (++$k>9) { $k = 0; shift(@unknown); }
      }
    }
  }


  sub create_position
  {
    my @data		= @_;
    my $p		= $data[1]." ".$data[2];
    my $k;

    for ($k=0; ($k<scalar(@data))&&(substr($data[$k],0,1) ne "#"); ++$k) { }
    @data		= splice(@data, $k-($k<8 ? 3 : 6), $k<8 ? 3 : 6);
    foreach (@L)
    {
      my @l		= split(" ");
      $p		= join(" ", $p, ($data[0]-$l[0])/$l[-1]+$data[3]);
      shift(@data);
    }
    return $p;
  }


  sub create_atom_ids
  {
    my $res		= 0;
    my $last		= 1;
    my @data;
    my $n;
    my $k;
    my $tmp;
    my $id;
    my %link;
    my %special;

    printf("Info: creating atom ids\n") if ($info);
    create_names();
    $n			= goto_data(atoms);
    foreach ((CL, NA, MG, K)) { $special{$_} = "ION"; }
    $special{"H H O"}	= "HOH";
    for (my $i=0; $i<=$n; ++$i)
    {
      @data		= split(" ", <DATA>) if ($i<$n);
      $traject[$i]	= create_position(@data) if ($i<$n);	# positions
      if ($i<$n ? ($mol[$i+1] = $data[1])!=$last : 1)		# residues
      {
	if ((($tmp = $link{$id = join(" ", sort(split(" ", $id)))}) eq "")&&
	    (($tmp = $special{$id}) eq ""))
	{
	  $link{$id}	=
	      $tmp 	= "R".($res<10 ? "0" : "").$res;
	  ++$res;
	}
	$cluster[$last]	= "PROT";
	$cluster[$last]	= "WATR" if ($tmp eq "HOH");
	$cluster[$last]	= "SALT" if ($tmp eq "ION");
	$residue[$last] = $tmp;
	$last		= $data[1];
	$id		= "";
      }
      $id		= join(" ", $id, $names[$data[2]]);
    }
  }


  sub crossover
  {
    my @d		= V_Subtr((split(" ", $position[@_[0]]))[0,1,2],
				  (split(" ", $position[@_[1]]))[0,1,2]);

    $d[0]		/= $l[3];
    $d[1]		/= $l[4];
    $d[2]		/= $l[5];
    return (($d[0]<-0.5)||($d[0]>=0.5)||($d[1]<-0.5)||($d[1]>=0.5)||
      ($d[2]<-0.5)||($d[2]>=0.5));
  }


  sub delete_crossovers
  {
    my $n		= scalar(@bonds);
    my $i		= 0;

    printf("Info: deleting crossover bonds\n") if ($info);
    while ($i<$n)					# take out crossovers
    {
      if (crossover(split(" ", $bonds[$i]))) { splice(@bonds, $i, 1); --$n; }
      else { ++$i; }
    }
  }


  sub delete_exclude
  {
    my $n		= scalar(@bonds);
    my $i		= 0;

    printf("Info: deleting excluded bonds\n") if ($info);
    while ($i<$n)
    {
      my $m		= $mol[(split(" ", $bonds[$i]))[0]+1];
      if (($exclude{$m} ne "")||($exclude{$residue[$m]} ne "")
	||($exclude{$cluster[$m]} ne "")) { splice(@bonds, $i, 1); --$n; }
      else { ++$i; }
    }
  }


  sub create_bonds
  {
    my $n		= goto_data(bonds);

    printf("Info: creating bonds\n") if ($info);
    for (my $i=0; $i<$n; ++$i)
    {
      my @arg		= split(" ", <DATA>);
      $bonds[$i]	= ($arg[2]-1)." ".($arg[3]-1);
    }
  }


# traject operators

  sub advance_traject
  {
    my $subject		= "item: ".lc(shift(@_));
    my $advance		= 1;

    while (!eof(TRAJECT)&&(substr(lc(join(" ", split(" ", <TRAJECT>))),
				  0,length($subject)) ne $subject)) {}
  }


  sub read_traject
  {
    my @box;
    my @l;

    advance_traject("timestep");
    my $timestep	= (split(" ", <TRAJECT>))[0];
    advance_traject("number of atoms");
    my $n		= (split(" ", <TRAJECT>))[0];
    advance_traject("box bounds");
    for (my $i=0; $i<3; ++$i)
    {
      my @data		= split(" ", <TRAJECT>);
      $box[$i]		= $data[0];			# box edge
      $l[$i]		= $data[1]-$data[0];		# box length
    }
    advance_traject("atoms");
    for (my $i=0; $i<$n; ++$i)
    {
      my @data		= split(" ", <TRAJECT>);	# read data
      $traject[$data[0]-1] = join(" ", @data);		# store data in order
    }
    return ($timestep, $n, @box, @l);
  }


# pdb format

  sub eigen_vector				# eigen_vector(@A, $l)
  {
    my @A		= splice(@_,0,9);
    my $l		= shift(@_);

    $A[0]		-= $l;
    $A[4]		-= $l;
    $A[8]		-= $l;
    @A			= M_Sweep(@A);
    return V_Norm(-$A[2], -$A[5], 1) if (($A[0]==1)&&($A[4]==1));
    return V_Norm(-$A[1], 1, -$A[7]) if (($A[0]==1)&&($A[8]==1));
    return V_Norm(1, -$A[3], -$A[6]) if (($A[4]==1)&&($A[8]==1));
    return (1,0,0) if ($A[0]==1);
    return (0,1,0) if ($A[4]==1);
    return (0,0,1) if ($A[8]==1);
    return (0,0,0);
  }


  sub pdb_inertia
  {
    my @s		= (
      @_[3]-@_[0]*@_[0], @_[4]-@_[1]*@_[1], @_[5]-@_[2]*@_[2],
      @_[6]-@_[0]*@_[1], @_[7]-@_[0]*@_[2], @_[8]-@_[1]*@_[2]);
    my @c		= (
      ($s[0]+$s[1])*(($s[0]+$s[2])*($s[1]+$s[2])-$s[3]*$s[3])		# c0
        -$s[5]*($s[3]*$s[4]+($s[1]+$s[2])*$s[5])
	-$s[4]*(($s[0]+$s[2])*$s[4]+$s[3]*$s[5]),
      -($s[0]+$s[2])*($s[1]+$s[2])-($s[0]+$s[1])*($s[0]+$s[1]+2*$s[2])	# c1
        +$s[3]*$s[3]+$s[4]*$s[4]+$s[5]*$s[5],
      2*($s[0]+$s[1]+$s[2]),						# c2
      -1);								# c3
    my @sol		= R_Third(@c);
    my @M		= ($s[1]+$s[2], -$s[3], -$s[4],
			   -$s[3], $s[0]+$s[2], -$s[5],
			   -$s[4], -$s[5], $s[0]+$s[1]);
    my @a		= eigen_vector(@M, $sol[0]);
    my @b		= eigen_vector(@M, $sol[2]);
    my @A		= M_Transpose(M_RotateNormal(@a));
    my @B		= M_Dot(M_RotateNormal(0,1,0),
			    M_Transpose(M_RotateNormal(MV_Dot(@A,@b))));
    return M_Dot(@B, @A);
  }



  sub pdb_focus		# using moment of inertia
  {
    my @R		= pdb_inertia(@_);

    printf("Info: focussing\n") if ($info);
    foreach (@position)
    {
      my @p		= split(" ");
      $_		= join(" ", MV_Dot(@R, @p[0,1,2]), @p[3,4]);
    }
  }


  sub pdb_center	
  {
    my @c		= splice(@_, 0, 3);

    printf("Info: centering\n") if ($info);
    @l[0,1,2]		= V_Mult(-1/2, @l[3,4,5]);
    foreach (@position)
    {
      my @p		= split(" ");
      $_		= join(" ", V_Subtr(@p[0,1,2], @c), @p[3,4]);
    }
  }


  sub pdb_pbc
  {
    printf("Info: applying periodicity\n") if ($info);
    foreach (@position)
    {
      my @p		= split(" ");
      my $m		= $mol[$p[3]];
      $_		= join(" ", V_PBC(@p[0,1,2], @l[3,4,5]), @p[3,4])
            if ($pbc{ALL}||$pbc{$m}||$pbc{$residue[$m]}||$pbc{$cluster[$m]});
    }
  }


  sub pdb_repair_bonds
  {
    return if (!$pbc);
    printf("Info: repairing bonds\n") if ($info);
    foreach (@bonds)
    {
      my @b		= split(" ");
      my @p		= split(" ", $position[$b[0]]);
      my @q		= split(" ", $position[$b[1]]);
      $position[$b[1]]	= join(" ", V_Add(@p[0,1,2], V_PBC(
			  V_Subtr(@q[0,1,2], @p[0,1,2]), @l[3,4,5])), @q[3,4]);
    }
  }


  sub pdb_positions
  {
    my @m		= (0,0,0,0,0,0,0,0,0);
    my $nmass		= 0;
    my $i		= 0;
    my $mass;
    my @p;
    my $d;

    foreach (@traject)
    {
      my @arg		= split(" ");
      my $m		= $mol[$arg[0]];
      next if (($exclude{$m} ne "")||($exclude{$residue[$m]} ne "")
	||($exclude{$cluster[$m]} ne ""));
      if ($units)
      {
	$p[0]		= $arg[2]+$arg[5]*$l[3];
	$p[1]		= $arg[3]+$arg[6]*$l[4];
	$p[2]		= $arg[4]+$arg[7]*$l[5];
      }
      else
      {
	$p[0]		= $l[0]+($arg[2]+$arg[5])*$l[3];
	$p[1]		= $l[1]+($arg[3]+$arg[6])*$l[4];
	$p[2]		= $l[2]+($arg[4]+$arg[7])*$l[5];
      }
      $position[$i++]	= join(" ", @p, $arg[0], $arg[1]);
      next if (($center{$m} eq "")&&($center{$cluster[$m]} eq ""));
      $nmass		+= $mass = $masses[$arg[1]];	# inertia necessities:
      $m[0]		+= $d = $mass*$p[0];		# <x>
      $m[3]		+= $d*$p[0];			# <x^2>
      $m[6]		+= $d*$p[1];			# <xy>
      $m[7]		+= $d*$p[2];			# <xz>
      $m[1]		+= $d = $mass*$p[1];		# <y>
      $m[4]		+= $d*$p[1];			# <y^2>
      $m[8]		+= $d*$p[2];			# <yz>
      $m[2]		+= $d = $mass*$p[2];		# <z>
      $m[5]		+= $d*$p[2];			# <z^2>
    }
    pdb_center(M_Mult(1/$nmass, @m)) if ($nmass);
    pdb_focus(M_Mult(1/$nmass, @m)) if ($nmass && scalar(%focus));
    pdb_pbc() if (scalar(%pbc));
    pdb_repair_bonds() if ($repair);
  }


  sub pdb_header
  {
    printf(PDB "REMARK  \n");
    printf(PDB "REMARK  ".$project."_trj.pdb GENERATED FROM ".$project.
      ($psf_ctrl ? "_ctrl.psf" : ".data")." AND $project.dump\n");
    printf(PDB "REMARK  CREATED BY $program v$version ON %s", `date`);
    printf(PDB "REMARK  \n");
    printf(PDB "CRYST1 ");
    printf(PDB "%8.8s ", $l[3]);
    printf(PDB "%8.8s ", $l[4]);
    printf(PDB "%8.8s ", $l[5]);
    printf(PDB "%6.6s ", "90");
    printf(PDB "%6.6s ", "90");
    printf(PDB "%6.6s ", "90");
    printf(PDB "%-11.11s", "P 1");
    printf(PDB "%3.3s\n", "1");
  }


  sub pdb_atoms
  {
    my $n		= 0;
    my $l		= 0;
    my $last		= 0;
    my @base;

    pdb_positions();
    psf_goto(atoms) if ($psf_ctrl);
    printf("Info: writing positions for timestep %d\n", $description[0])
      if ($info);
    foreach (@position)
    {
      my @p		= split(" ");
      my $nres		= $mol[$p[3]];
      my @psf		= split(" ", advance_read(PSF_CTRL, $p[3]-$l))
							      if ($psf_ctrl);
      @base		= @p[0,1,2] if ($last!=$nres);
      #@p		= V_Add(V_PBC(V_Subtr(@p, @base), @l[3,4,5]), @l);
      foreach (@p) { $_ = 0 if (abs($_)<1e-4); }
      printf(PDB "ATOM  ");				# pdb command
      printf(PDB "%6.6s ", ++$n);			# atom number
      printf(PDB "%-3.3s ",
	$psf_ctrl ? $psf[4] : $names[$p[4]]);		# atom name
      printf(PDB "%-3.3s ",
	$psf_ctrl ? $psf[3] : $residue[$nres]);		# residue name
      printf(PDB "%5.5s ", $nres);			# residue number
      printf(PDB "%3.3s ", "");				# empty placeholder
      printf(PDB "%7.7s %7.7s %7.7s ",
       	$p[0], $p[1], $p[2]);				# positions
      printf(PDB "%5.5s %5.5s %4.4s ",
	"1.00", "0.00", "");				# trailing variables
      printf(PDB "%-4.4s\n",
       	$psf_ctrl ? $psf[1] : $cluster[$nres]);		# cluster name
      $last		= $nres;
      $l		= $p[3];
    };
    printf(PDB "END\n");
  }


  sub pdb_timestep
  {
    pdb_header();
    pdb_atoms();
    return ();
  }


# psf format

  sub psf_header
  {
    printf(PSF "PSF\n");
    printf(PSF "\n");
    printf(PSF "%8.8s !NTITLE\n", 2);
    printf(PSF "REMARK   ".$project."_trj.psf GENERATED FROM $project.data\n");
    printf(PSF "REMARK   CREATED BY $program v$version ON %s", `date`);
    printf(PSF "\n");
  }


  sub psf_atoms
  {
    my $n		= goto_data(atoms);
    my $natom		= 0;
    my $l		= 0;
    my $k		= 0;
    my @extra;

    for (my $i=0; $i<$n; ++$i)
    {
      my @arg		= split(" ", <DATA>);
      if (!$k) {
	for ($k=0; ($k<scalar(@arg))&&(substr($arg[$k],0,1) ne "#"); ++$k) {} }
      $extra[$arg[0]-1]	= $arg[1]." ".$arg[2]." ".$arg[3];
    }
    printf(PSF "%8.8s !NATOM\n", scalar(@position));
    foreach (@position)
    {
      my @p		= split(" ");
      my @arg		= split(" ", $extra[$natom]);
      printf(PSF "%8.8s ", ++$natom);			# atom number
      printf(PSF "%-4.4s ", $cluster[$arg[0]]);		# cluster name
      printf(PSF "%-4.4s ", $arg[0]);			# residue number
      printf(PSF "%-4.4s ", $residue[$arg[0]]);		# residue name
      printf(PSF "%-4.4s ", $names[$p[4]]);		# atom name
      printf(PSF "%4.4s ", $arg[1]);			# atom number
      printf(PSF "%10.10s ", $k%3 ? $arg[2] : 0);	# atom charge
      printf(PSF "%4.4s ", "");				# blank entry
      printf(PSF "%8.8s ", $masses[$arg[1]]);		# atom mass
      printf(PSF "%11.11s\n", "0");			# trailing variable
      $l		= $p[3];
      last if ($natom>=$n)
    }
    printf(PSF "\n");
  }


  sub psf_bonds
  {
    my $npairs		= 0;

    delete_exclude() if (scalar(%exclude)>0);
    delete_crossovers() if ($cut);
    printf(PSF "%8.8s !NBOND\n", scalar(@bonds));
    foreach (@bonds)
    {
      my @arg		= split(" ");
      printf(PSF " %7.7s %7.7s", $arg[0]+1, $arg[1]+1);
      if (++$npairs>=4) { printf(PSF "\n"); $npairs = 0; }
    }
    if ($npairs>0) { printf(PSF "\n"); }
    printf(PSF "\n");
  }


# main

  initialize();

  # create .pdb file

  $ncurrent		= -1;
  while ($traject_flag&&!eof(TRAJECT))
  {
    @description	= read_traject();
    @l			= splice(@description, -6, 6);
    next if (($nstart<0)||(++$ncurrent<$nstart));
    if ($dn<1) { pdb_timestep(); last; }
    $ncurrent		= pdb_timestep() if ($nstart||!($ncurrent%$dn));
    $nstart		= 0;
  }
  pdb_timestep() if ($nstart<0);

  # create .psf file

  if (!$psf_ctrl)
  {
    psf_header();
    psf_atoms();
    psf_bonds();
  }
  else
  {
    system("rm -f ".$psf_name);
    system("ln -s ".$psf_ctrl_name." ".$psf_name);
  }

  # add tail to files

  #printf(PDB "END\n");
  printf("\n") if ($info);
  close(PDB);
  close(PSF);

  close(TRAJECT);
  close(DATA);

