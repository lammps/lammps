#!/usr/bin/perl
# Computes potential energy of atom as a function of distance from another atom 
# and computes numerical derivates of potential. 
# The script was used to check if results from LAMMPS (using 2particles.in)
# are the same as these computed b this script.
# Prints results to STDOUT.
# Hence, use it like this:
# ./LebDer.pl > PerlResult.dat
# After that use lebedeva00.plot
#
# Author: Zbigniew Koziol, National Center for Nuclear Research, Poland
# Email: softquake@gmail.com

# Parameters used by ZJK for Lebedeva
my $LEB_A = -14.558;
my $LEB_B = 21.204;
my $LEB_alpha = 4.16;
my $LEB_C = 1.8;
my $LEB_D1 = -0.862;
my $LEB_D2 = 0.10049; # has very strong influence on position of minimum
my $LEB_lambda1 = 0.6; # has influance on splitting of AB-AA.
my $LEB_lambda2 = 0.4; # has strong influence on position of minimum
my $LEB_z0 = 3.198;
my $LEBSCALE =1.0;

$Z0=3.35;

$CX0 = 10;
$CY0 = 10;

for (my $t=0; $t<400; $t++) { 
	my $X0 = 0.001 + 0.05*$t;
	my $Y0 = 0.001 + 0.05*$t;
	my $Z  = $Z0;
	print $X0, "\t", $Y0, "\t", $Z, "\t",&LEB($X0, $Y0, $Z), "\t", &DLEBX($X0, $Y0, $Z),"\t",&DLEBY($X0, $Y0, $Z), "\t", &DLEBZ($X0, $Y0, $Z),"\n";
}

###############################################################################################

sub LEB {
	my $x = shift;
	my $y = shift;
	my $z = shift;
	
	my $rho2 = ($x-$CX0)*($x-$CX0) + ($y-$CY0)*($y-$CY0);
	my $r   = sqrt($rho2 + ($Z0)*($Z0));
	my $zr  = ($LEB_z0/$r)*($LEB_z0/$r);
	my $zr6 = $zr*$zr*$zr;

	my $ONE = $LEB_C*(1+$LEB_D1*$rho2+$LEB_D2*$rho2*$rho2);
	my $TWO = exp(-$LEB_lambda1*$rho2)*exp(-$LEB_lambda2*($z*$z-$LEB_z0*$LEB_z0));
	my $U   = $LEB_A*$zr6 +$LEB_B*exp(-$LEB_alpha*($r-$LEB_z0)) + $ONE*$TWO;
	return $U;
}

sub DLEBX { # finding derivative at $x
	my $x = shift;
	my $y = shift;
	my $z = shift;

        my $h = 0.0001;

        my $D = (&LEB($x+$h, $y, $z)-&LEB($x-$h, $y, $z))/(2*$h);
        
        return $D;
}

sub DLEBY { # finding derivative at $y
	my $x = shift;
	my $y = shift;
	my $z = shift;

        my $h = 0.0001;

        my $D = (&LEB($x, $y+$h, $z)-&LEB($x, $y-$h, $z))/(2*$h);
        
        return $D;
}

sub DLEBZ { # finding derivative at $z
	my $x = shift;
	my $y = shift;
	my $z = shift;

        my $h = 0.0001;

        my $D = (&LEB($x, $y, $z+$h)-&LEB($x, $y, $z-$h))/(2*$h);
        
        return $D;
}
