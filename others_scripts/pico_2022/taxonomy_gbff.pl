#!/usr/local/bin/perl
use File::Copy;



while ($file = glob("./gbff_files/*.gbff")) {

##############################################################
###                                                        ###
###  AQUI VAMOS A EXTRAER EL NOMBRE DE CADA ARCHIVO GBFF   ###
###                                                        ###
##############################################################

	open (LAL, ">>GBFF.NAMES.tmp") or die ("No puedo abrir GBFF.NAMES.tmp");
		print LAL ("$file\n");
	close (LAL);
}

open (LIZ, "GBFF.NAMES.tmp") or die ("No puedo abrir GBFF.NAMES.tmp");
	while (my $file = <LIZ>){
		$spp = $file;
		$spp =~ s/genomic\.gbff//g; #./gbff_files/Cyanobium_gracile_PCC_6307_PCC_6307__genomic.gbff
		$spp =~ s/\.\/gbff_files\///g;
		$spp =~ s/\n//g;
		open (LAL, $file) or die ("No puedo abrir $file\n");
			open (KAR, ">>taxonomy.names.txt") or die ("No puedoabrir species.names.txt\n");
				while (my $line = <LAL>){
					chomp ($line);
						if ($line =~/(Bacteria; .+;\n.+)/){# Bacteria; .+;\n.+
						$ORGANISM = $1;
						}
						if ($line =~ /.+(Bacteria); (.+); (.+); (.+);/){
								$name1=$1;
								$name1 =~ s/\n//g;
								$name2=$2;
								$name2 =~ s/\n//g;
								$name3=$3;
								$name3 =~ s/\n//g;
								$name4=$4;
								$name4 =~ s/\n//g;
								$name5=$5;
								$name5 =~ s/\n//g;
								#print("$spp\n");
								#print ("$name1\n");
								#print ("$name2\n");
								#print ("$name3\n");
								#print ("$name4\n");
								#print ("$name5\n");
						}
					}
				print KAR ("$spp\t$name4\n");
				print("$ORGANISM\n");
			close (KAR);
		close (LAL);
	}
close (LIZ);

unlink glob("./*.tmp");
#unlink glob("./*.txt");

