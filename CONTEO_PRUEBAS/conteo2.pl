#!/usr/local/bin/perl

use strict;
my $count_pal =0;
my $count_A = 0;
my $count_C = 0;
my $count_G = 0;
my $count_T = 0;
my $count_GC = 0;
my $count_CG = 0;
my $conteo = 0;

open (LUZ, "sequence.txt") or die ("sequence.txt\n");
	while( <LUZ> ) {
	   #print;
	    ++$count_pal while m[GCGC]ig;
	    ++$count_A while m[A]ig;
	    ++$count_C while m[C]ig;
	    ++$count_G while m[G]ig;
	    ++$count_T while m[T]ig;
	    ++$count_GC while m[GC]ig;
	    ++$count_CG while m[CG]ig;
	}
	while (my $line = <LUZ>){
		chomp ($line);
		$conteo = $conteo+(length($line));
				}
	my $EXP = (($count_GC*$count_CG*$count_GC)/($count_C*$count_G))*($conteo-4+1);
	
	print "'A' appeared $count_A times in 'the File'\n";
	print "'C' appeared $count_C times in 'the File'\n";
	print "'G' appeared $count_G times in 'the File'\n";
	print "'T' appeared $count_T times in 'the File'\n";
	print "EXP $EXP\n";
	print "OBS $count_pal\n";
	print "LONGITUD $conteo\n";
close (LUZ);
