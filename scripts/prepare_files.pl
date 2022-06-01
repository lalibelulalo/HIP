#!/usr/local/bin/perl

my $o;
my $l = 0;
my $name;

print ("				*	*	*\n\nEste script descarga los archivos fasta (faa & fna) para cada accesion en el genoma.\n\n				*	*	*\n\n");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

my $lastdate = "$mday-$mon-$year\_$hour\:$min";
system("mkdir Archivos_$lastdate/
	mkdir Genome_fasta_aa/
	mkdir Genome_fasta_nu/
	mkdir ACCESIONES");
	
=pod
open (PAM, ">FALTAN.faa.txt") or die ("No puedo abrir FALTAN.faa.txt\n");
	print ALE ("Especie\tAccesion\n");
close (PAM);
open (PAM, ">FALTAN.fna.txt") or die ("No puedo abrir FALTAN.fna.txt\n");
	print ALE ("Especie\tAccesion\n");
close (PAM);
=cut
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
my $t=0;
open (MAR, "species.names") or die ("No puedoabrir species.names\n");
	while (my $line = <MAR>){
		chomp ($line);
		print ("Especie ",$t+1," : ",$line,"\n");
		##############################################
		##                                          ##
		##          OBTENIENDO ACCESIONES           ##
		##                                          ##
		##############################################
		my $l = 0;
		my $z = 0;

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
		system ("mkdir Genome_fasta_aa_$line/
			 mkdir Genome_fasta_nu_$line");
		#system ("mkdir Genbank_aa");
		
		open (LIZ, "$line.ACCESIONES.txt") or die ("No puedo abrir $line.ACCESIONES.txt\n");
		while (my $linea = <LIZ>){
			chomp ($linea);
			if ($linea =~ /(\w+)/){
				my $acc = $1;
				print ("	Accesion ",$z+1," : ", "$acc\n");
				system ("esearch -db nuccore -query $acc | elink -target protein | efetch -format fasta >$acc.faa");
				system ("esearch -db nuccore -query $acc | efetch -format fasta >$acc.fna");
				
				if (-z "$acc.faa") {
					print ("		Falta: $line \t$acc.faa\n");
					open (PAM,'>>',"$line.FALTAN.faa.txt") or die ("No puedo abrir $line.FALTAN.faa.txt\n");
						print PAM ("$line\t$acc\n");
					close (PAM);
						}
				
				if (-z "$acc.fna") {
					print ("		Falta: $line \t$acc.fna\n");
					open (ALE,'>>',"$line.FALTAN.fna.txt") or die ("No puedo abrir $line.FALTAN.fna.txt\n");
						print ALE ("$line\t$acc\n");
					close (ALE);
						}
				system("mv *.faa Genome_fasta_aa_$line");
				system("mv *.fna Genome_fasta_nu_$line");
						
				if (-e "$line.FALTAN.faa.txt") {
								system("cat *.faa.txt >>FALTAN.aa.txt");
								system("rm $line.FALTAN.faa.txt");
				}
				if (-e "$line.FALTAN.fna.txt") {
								system("cat *.fna.txt >>FALTAN.nu.txt");
								system("rm $line.FALTAN.fna.txt");				
				}
				
				sleep(20);
				$z++;
				
			}
			


		}
		close (LIZ);	
		$t++;
		#system("cat *.faa.txt >FALTAN.faa.txt");
		#system("cat *.fna.txt >FALTAN.faa.txt");			
				
		system("mv Genome_fasta_nu_$line Genome_fasta_nu/
		mv Genome_fasta_aa_$line Genome_fasta_aa/
		mv $line.ACCESIONES.txt ACCESIONES");
		
		
		
				
	}

close (MAR);

if (-e "FALTAN.aa.txt") {
				system("mv FALTAN.aa.txt Archivos_$lastdate");
			}
		
if (-e "FALTAN.nu.txt") {
				system("mv FALTAN.nu.txt Archivos_$lastdate");	
			}

system("mv Genome_fasta_aa Archivos_$lastdate/
	mv Genome_fasta_nu Archivos_$lastdate/
	mv ACCESIONES Archivos_$lastdate/
	mv FASTA.NAMES Archivos_$lastdate/
	mv species.names Archivos_$lastdate");

