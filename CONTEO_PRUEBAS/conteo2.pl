#!/usr/local/bin/perl
my $o;
my $l = 0;
my $name;

############################################################################################
########          ESTE ARCHIVO VA A IR GUARDANDO LOS RESULTADOS DEL CONTEO         #########
############################################################################################
open (MAR, ">RESULTADOS_ABUNDANCIA.csv") or die ("No puedo abrir RESULTADOS_ABUNDANCIA.csv\n");
	print MAR ("species,accesion,GC,CG,GA,AT,TC,A,T,G,C,LEN,EXP,OBS\n");
close (MAR);


open (MAR, "species.names") or die ("No puedoabrir species.names\n"); #VAMOS IR ENTRANDO A CADA CARPETA DE CADA ESPECIE
	while (my $line = <MAR>){		
		chomp ($line);
		print ("Especie ",$o+1,": ",$line,"\n");
		$o++;
		my $z = 0;
		open (LIZ, "ACCESIONES/$line.ACCESIONES.txt") or die ("No puedo abrir $line.ACCESIONES.txt\n"); #EN CADA CARPETA DE ESPECIE VAMOS IR LEYENDO CADA ARCHIVO
		while (my $linea = <LIZ>){
			chomp ($linea);
			if ($linea =~ /(\w+)/){
				my $acc = $1;
				print ("\tAccesion ",$z+1,": ",$acc,"\n");
				
				system ("sed '/>/d' Genome_fasta_nu/Genome_fasta_nu_$line/$acc.fna >$line.$acc.fna.headless"); #LE QUITAMOS EL ENCABEZADO PARA PODER HACER EL CONTEO UNICAMENTE CON LOS NUCLEOTIDOS
				
				open (GEM, ">>RESULTADOS_ABUNDANCIA.csv") or die ("No puedo abrir RESULTADOS_ABUNDANCIA.csv\n");	
					my $conteo;
					open (LUZ, "$line.$acc.fna.headless") or die ("$line.$acc.fna.headless\n"); #EN ESTA LINEA HACEMOS EL CONTEO DE TODOS LOS NUCLEOTIDOS
						while (my $line = <LUZ>){
							chomp ($line);
							$conteo = $conteo+(length($line));		
						}
					close (LUZ);
					#AQUI VAMOS A HACER EL CONTEO POR NUCLEOTIDO, DINUCLEOTIDO.		
					$gc = `grep -c 'GC' $line.$acc.fna.headless`;
					$gc =~ s/\n//g;
					$cg = `grep -c 'CG' $line.$acc.fna.headless`;
					$cg =~ s/\n//g;
					$ga = `grep -c 'GA' $line.$acc.fna.headless`;
					$ga =~ s/\n//g;
					$at = `grep -c 'AT' $line.$acc.fna.headless`;
					$at =~ s/\n//g;
					$tc = `grep -c 'TC' $line.$acc.fna.headless`;
					$tc =~ s/\n//g;
					$a = `grep -c 'A' $line.$acc.fna.headless`;
					$a =~ s/\n//g;
					$t = `grep -c 'T' $line.$acc.fna.headless`;
					$t =~ s/\n//g;
					$g = `grep -c 'G' $line.$acc.fna.headless`;
					$g =~ s/\n//g;
					$c = `grep -c 'C' $line.$acc.fna.headless`;
					$c =~ s/\n//g;
					$pal = `grep -c 'GCGATCGC' $line.$acc.fna.headless`;#ABUNDANCIA HIP1
					$EXP = (($gc*$cg*$ga*$at*$tc*$cg*$gc)/($c*$g*$a*$t*$c*$g))*($conteo-8+1);#ABUNDANCIA ESPERADA (EXP)
					my $renglon = "$line,$linea,$gc,$cg,$ga,$at,$tc,$a,$t,$g,$c,$conteo,$EXP,$pal";
					
					print GEM ("$renglon");
					system ("rm $line.$acc.fna.headless");
				close (GEM);	
				
				sleep(20);
				$z++;
			}
		}	
		close (LIZ);
	}
close (MAR);
