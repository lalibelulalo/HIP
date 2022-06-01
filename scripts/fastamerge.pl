#!/usr/bin/perl -w
my $o;

system("mkdir all_merged");

open (MAR, "species.names") or die ("No puedoabrir species.names\n");
	while (my $line = <MAR>){
		chomp ($line);
		$o++;
		my $l = 0;
		my $name;
		
		print ("Carpeta: Genome_fasta_aa_$line\n");
		system("cat Genome_fasta_aa_$line/*.faa > Genome_fasta_aa_$line/$line.merged.faa");
		system("mv Genome_fasta_aa_$line/$line.merged.faa all_merged");
		print ("Siguente archivo\n");

	}
close (MAR);

system("orthofinder -f all_merged");
