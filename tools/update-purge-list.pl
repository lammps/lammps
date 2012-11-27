#!/usr/bin/perl -w

use strict;
my $date="";
my $newdate="";
my @list;
my $scm=0;

open(GITLOG,'git log --summary |');
open(LIST,'> purge-list.txt');

while (<GITLOG>) {
  chomp;
  if (/^Date:\s+(.*)/) {
    $newdate = $1;
  }

  if (/\s+delete mode \d+ src\/(\S+)\/(\S+\.(cpp|h)).*/) {
    # check if file exists in a different sub directory
    @list = glob("[A-Z-][A-Z-]*/$2");
    if ($#list < 0) {
      # check if file got moved to main source directory
      $scm = system("git ls-files $2 --error-unmatch < /dev/null > /dev/null 2>&1");
      if ($scm) {
        if ($date ne $newdate) {
          $date = $newdate;
          print LIST "# deleted on $date\n";
        }
        print LIST "$2\n";
      }
    }
  }
}

close GITLOG;
close LIST;
