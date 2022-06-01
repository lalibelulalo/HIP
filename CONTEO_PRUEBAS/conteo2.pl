#!/usr/local/bin/perl

use strict;

my $count = 0;
open (LUZ, "sequence.txt") or die ("sequence.txt\n");
	while( <LUZ> ) {
	   #print;
	    ++$count while m[A]ig;
	}
	print "'A' appeared $count times in 'the File'\n";
close (LUZ);
