#!/usr/local/bin/perl

my $conteo;
open (MAR, "NZ_CP011382.fna.headless") or die ("No puedoabrir NZ_CP011382.fna.headless\n");
	while (my $line = <MAR>){
		chomp ($line);
		$conteo = $conteo+(length($line));		
	}
	print ("El conteo final es: ", $conteo,"\n");
close (MAR);
