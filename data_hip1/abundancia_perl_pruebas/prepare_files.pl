#!/usr/local/bin/perl

my $o;
my $l = 0;
my $name;

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

my $lastdate = "$mday-$mon-$year\_$hour\:$min";
system("mkdir Archivos_$lastdate");

while ($file = glob("./*.fasta")) {

##############################################################
###                                                        ###
###  AQUI VAMOS A EXTRAER EL NOMBRE DE CADA ARCHIVO FASTA  ###
###                                                        ###
##############################################################

	open (LAL, ">>FASTA.NAMES") or die ("No puedo abrir FASTA.NAMES");
		print LAL ("$file\n");
	close (LAL);
}


open (LAL, "FASTA.NAMES") or die ("No puedo abrir FASTA.NAMES\n");
	open (KAR, ">species.names") or die ("No puedoabrir species.names\n");
		while (my $line = <LAL>){
			chomp ($line);
			$o++;
			$line =~ s/\.\///g;
			$line =~ s/\.fasta//g;
		print KAR ("$line\n");
		}
	close (KAR);
close (LAL);

##############################################################
###                                                        ###
###  	    OBTENIENDO ACCESIONES DE CADA GENOMA	   ###
###                                                        ###
##############################################################
open (MAR, "species.names") or die ("No puedoabrir species.names\n");
	while (my $line = <MAR>){
		chomp ($line);
		$o++;
		##############################################
		##                                          ##
		##          OBTENIENDO ACCESIONES           ##
		##                                          ##
		##############################################
		my $l = 0;
		my $name;

		open (AUR, "$line.fasta") or die ("No puedo abrir $line.fasta\n");
			open (LIZ, ">$line.ACCESIONES.txt") or die ("No puedo abrir $line.ACCESIONES.txt\n");
			while (my $linea = <AUR>){
				chomp ($linea);
				$l++;
				
				if ($linea =~ />(\w+) /){#>NC_019727
					$name=$1;
					print LIZ ("$name\n");
				} 
					
			}
			close (LIZ);
		close (AUR);
			
		##########################################################################
		###                                                                    ###
		###  Obtenemos los archivos Genbank de las proteínas correspondientes  ###
		###                                                                    ###
		##########################################################################
		
		#system ("mkdir Genbank_aa_$line");
		system ("mkdir Genome_fasta_aa_$line");
		system ("mkdir Genome_fasta_nu_$line");
		#system ("mkdir Genbank_aa");
		system ("mkdir Genome_fasta_aa");
		system ("mkdir Genome_fasta_nu");
		system ("mkdir ACCESIONES");
		
		my $z = 0;

		open (LIZ, "$line.ACCESIONES.txt") or die ("No puedo abrir $line.ACCESIONES.txt\n");
		while (my $linea = <LIZ>){
			chomp ($linea);
			if ($linea =~ /(\w+)/){
				my $acc = $1;
				print ("$acc\n");
				#my $outfile1 = $acc.'.aa.gbff';
				my $outfile2 = $acc.'.faa';
				my $outfile3 = $acc.'.fna';
				#system ("curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=gb&retmode=text&id=$acc' > Genbank_aa_$line/$outfile1");
				system ("esearch -db nuccore -query $acc | elink -target protein | efetch -format fasta > Genome_fasta_aa_$line/$outfile2");
				system ("esearch -db nuccore -query $acc | efetch -format fasta > Genome_fasta_nu_$line/$outfile3");

				sleep(20);
				$z++;
				print ("archivo numero: $z\n");
			}
			system("mv Genome_fasta_aa_$line Genome_fasta_aa");
			system("mv Genome_fasta_nu_$line Genome_fasta_nu");
			system("mv $line.ACCESIONES.txt ACCESIONES");
			
		}
		close (LIZ);	
		
		#system("esearch -db nuccore -query NZ_AP018172 | elink -target protein | efetch -format fasta >NZ_AP018172.faa");
		
	}

close (MAR);

system("mv Genome_fasta_aa Archivos_$lastdate");
system("mv Genome_fasta_nu Archivos_$lastdate");
system("mv ACCESIONES Archivos_$lastdate");
system("mv FASTA.NAMES Archivos_$lastdate");
system("mv species.names Archivos_$lastdate");

