#!/usr/local/bin/perl

use strict;
=pod
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
=cut
my $o = 0;
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
	print GEM ("ARCH,A,f_A,T,f_T,C,f_C,G,f_G,GC,f_GC,CG,f_CG,GA,f_GA,AT,f_AT,TC,f_TC,LONG,f_obs,f_exp,num_pal_obs\n");
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
		# EN CADA CARPETA DE ESPECIE VAMOS IR LEYENDO CADA ARCHIVO
		open (LIZ, "ACCESIONES/$especie.ACCESIONES.txt") or die ("No puedo abrir $especie.ACCESIONES.txt\n"); 
		while (my $accesion = <LIZ>){
			chomp ($accesion);
			if ($accesion =~ /(\w+)/){
				# CONTADORES EN 0'S PARA EL CONTEO EN LA SIGUIENTE ESPECIE
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

				my $acc = $1;
				print ("\tAccesion ",$z+1,": ",$acc,"\n");
				# LE QUITAMOS EL ENCABEZADO PARA PODER HACER EL CONTEO UNICAMENTE CON LOS NUCLEOTIDOS
				system ("sed '/>/d' Genome_fasta_nu/Genome_fasta_nu_$especie/$acc.fna >$especie.$acc.fna.headless");
				
				# QUITAMOS LOS SALTOS DE LINEA PARA EVITAR SESGOS EN EL CONTEO
				system ("tr -d '\n' <$especie.$acc.fna.headless >$especie.$acc.fna.headless.1row");




				# AQUI SE HACE EL CONTEO
				#---------------------------------------------------------------------------
						# LONGITUD DE SECUENCIA
						open (BEL, "$especie.$acc.fna.headless.1row") or die ("No puedo abrir $especie.$acc.fna.headless.1row\n");
							while (my $line = <BEL>){
									chomp ($line);
									$long_genoma = $long_genoma+(length($line));			
									}
						close (BEL);


						open (CEC, "$especie.$acc.fna.headless.1row") or die ("No puedo abrir $especie.$acc.fna.headless.1row\n");
							while( <CEC> ) {
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
							# FRECUENCIA RELATIVA DE NUCLEOTIDOS Y DINUCLEOTIDOS
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
							#print "'T' aparece $count_T veces en el archivo '$especie.$acc.fna.headless'\n";
							#print "'C' aparece $count_C veces en el archivo '$especie.$acc.fna.headless'\n";
							#print "'G' aparece $count_G veces en el archivo '$especie.$acc.fna.headless'\n\n";
							print "OBS: $count_pal\n";
							print "EXP: $EXP\n";
							print "LONGITUD $long_genoma\n";
							print "-----------------------------------------------\n";
							
							# GUARDANDO CONTEO EN RES.CONTEO.csv			
							open (MIR,">>RES.CONTEO.csv") or die ("No puedo abrir RES.CONTEO.csv");
								print MIR ("$especie.$acc,$count_A,$A,$count_T,$T,$count_C,$C,$count_G,$G,$count_GC,$GC,$count_CG,$CG,$count_GA,$GA,$count_AT,$AT,$count_TC,$TC,$long_genoma,$PAL,$EXP,$count_pal\n");
							close (MIR);
						close (CEC);
				#---------------------------------------------------------------------------	
				#AQUI TERMINA EL CONTEO
				system ("rm $especie.$acc.fna.headless");
				system ("rm $especie.$acc.fna.headless.1row");
				
				$z++;
			}
		}	
		close (LIZ);
	}
close (MAR);
