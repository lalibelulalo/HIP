#!/usr/bin/perl -w
my $o;


open (MAR, "species.names") or die ("No puedoabrir species.names\n");
	while (my $line = <MAR>){
		chomp ($line);
		$o++;
		my $l = 0;
		my $name;
		
		print ("Carpeta: Genome_fasta_aa_$line\n");
		system(" Genome_fasta_aa_$line/$line.faa");
		print ("Siguente archivo:\n");

	}

close (MAR);

