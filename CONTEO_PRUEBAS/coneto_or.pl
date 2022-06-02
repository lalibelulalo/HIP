#!/usr/local/bin/perl

use strict;
my $o = 0;
my $long_genoma = 0;
my $count_pal = 0;
my $count_A = 0;
my $count_T = 0;
my $count_C = 0;
my $count_G = 0;
my $count_GC = 0;
my $count_CG = 0;
my $count_GA = 0;
my $count_AT = 0;
my $count_TC = 0;
my $PAL = 0;
my $A = 0;
my $T = 0;
my $C = 0;
my $G = 0;
my $AA = 0;
my $AT = 0;
my $TA = 0;
my $TT = 0;
my $GC = 0;
my $CG = 0;
my $GA = 0;
my $AT = 0;
my $TC = 0;

# CREANDO ARCHIVO DE RESULTADOS
open (GEM,">RES.CONTEO.csv") or die ("No puedo abrir RES.CONTEO.csv");
	print GEM ("ARCH,f_A,f_T,f_C,f_G,f_GC,f_CG,f_GA,f_AT,f_TC,LONG,f_obs,f_exp,num_pal_obs\n");
			#$A,$T,$C,$G,$GC,$CG,$GA,$AT,$TC
close (GEM);



# EXTRAYENDO NOMBRES DE ARCHIVOS
while (my $file = glob("./*_600.txt")) {

	open (LAL, ">>SEQ.NAMES") or die ("No puedo abrir SEQ.NAMES");
		$file =~ s/\.\///g;
		print LAL ("$file\n");
	close (LAL);
}



open (MAR, "species.names") or die ("No puedoabrir species.names\n"); #VAMOS IR ENTRANDO A CADA CARPETA DE CADA ESPECIE
	while (my $especie = <MAR>){		
		chomp ($especie);
		print ("Especie ",$o+1,": ",$especie,"\n");
		$o++;
		my $z = 0;
		#EN CADA CARPETA DE ESPECIE VAMOS IR LEYENDO CADA ARCHIVO
		open (LIZ, "ACCESIONES/$especie.ACCESIONES.txt") or die ("No puedo abrir $especie.ACCESIONES.txt\n"); 
		while (my $accesion = <LIZ>){
			chomp ($accesion);
			if ($accesion =~ /(\w+)/){
				my $acc = $1;
				print ("\tAccesion ",$z+1,": ",$acc,"\n");
				#LE QUITAMOS EL ENCABEZADO PARA PODER HACER EL CONTEO UNICAMENTE CON LOS NUCLEOTIDOS
				system ("sed '/>/d' Genome_fasta_nu/Genome_fasta_nu_$especie/$acc.fna >$especie.$acc.fna.headless");
				#system ("z 's/\n/,/g;s/,$/\n/' $especie.$acc.fna.headless >$especie.$acc.1row.fna.headless");




				# AQUI SE HACE EL CONTEO
				#---------------------------------------------------------------------------
						# LONGITUD DE SECUENCIA
						open (BEL, "$especie.$acc.1row.fna.headless") or die ("No puedo abrir $especie.$acc.1row.fna.headless\n");
							while (my $line = <BEL>){
									chomp ($line);
									$long_genoma = $long_genoma+(length($line));
									#print ("El conteo va: $long_genoma\n");			
									}
						close (BEL);


						open (CEC, "$especie.$acc.1row.fna.headless") or die ("No puedo abrir $especie.$acc.1row.fna.headless\n");
							while( <CEC> ) {
							#print;
								++$count_pal while m[GCGATCGC]ig;
								++$count_A while m[A]ig;
								++$count_T while m[T]ig;
								++$count_C while m[C]ig;
								++$count_G while m[G]ig;
								++$count_GC while m[GC]ig;
								++$count_CG while m[CG]ig;
								++$count_GA while m[GA]ig;
								++$count_AT while m[AT]ig;
								++$count_TC while m[TC]ig;
							}
							$PAL = $count_pal/($long_genoma-4+1);
							$A = $count_A/($long_genoma-1+1);
							$T = $count_T/($long_genoma-1+1);
							$C = $count_C/($long_genoma-1+1);
							$G = $count_G/($long_genoma-1+1);
							$GC = $count_GC/($long_genoma-2+1);
							$CG = $count_CG/($long_genoma-2+1);
							$GA = $count_GA/($long_genoma-2+1);
							$AT = $count_AT/($long_genoma-2+1);
							$TC = $count_TC/($long_genoma-2+1);
							
							my $EXP = (($GC*$CG*$AT*$TC*$CG*$GC)/($C*$G*$A*$T*$C*$G))*($long_genoma-8+1);
							my $count_pal_EXP = (($GC*$CG*$AT*$TC*$CG*$GC)/($C*$G*$A*$T*$C*$G));
							print "-----------------------------------------------\n";
							#print "'A' aparece $count_A veces en el archivo '$especie.$acc.fna.headless'\n";
							#print "'T' aparece $count_T veces en el archivo '$especie.$acc.fna.headless'\n\n";
							#print "EXP $EXP\n";
							print "OBS: $count_pal\n";
							print "EXP: $EXP\n";
							print "LONGITUD $long_genoma\n";
							print "-----------------------------------------------\n";
										
							open (MIR,">>RES.CONTEO.csv") or die ("No puedo abrir RES.CONTEO.csv");
								print MIR ("$especie.$acc,$A,$T,$C,$G,$GC,$CG,$GA,$AT,$TC,$long_genoma,$PAL,$EXP,$count_pal\n");
							close (MIR);
							
						my $count_pal = 0;
						my $count_A = 0;
						my $count_T = 0;
						my $count_C = 0;
						my $count_G = 0;
						my $count_GC = 0;
						my $count_CG = 0;
						my $count_GA = 0;
						my $count_AT = 0;
						my $count_TC = 0;
						my $PAL = 0;

						close (CEC);
						$long_genoma = 0;
				#---------------------------------------------------------------------------	
				#AQUI TERMINA EL CONTEO
				system ("rm $especie.$acc.1row.fna.headless");
				
				$z++;
			}
		}	
		close (LIZ);
	}
close (MAR);

