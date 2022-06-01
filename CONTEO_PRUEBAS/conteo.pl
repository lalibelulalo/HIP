#!/usr/local/bin/perl

	
#AQUI VAMOS A HACER EL CONTEO POR NUCLEOTIDO, DINUCLEOTIDO.
my $conteo=0;
open (LUZ, "sequence.txt") or die ("sequence.txt"); #EN ESTA LINEA HACEMOS EL CONTEO DE TODOS LOS NUCLEOTIDOS
						while (my $line = <LUZ>){
							chomp ($line);
							$conteo = $conteo+(length($line));
							}	

		
$gc = system("grep -o 'GC' sequence.txt | wc -l");
$gc =~ s/\n//g;
$cg = system("grep -o 'CG' sequence.txt | wc -l");
$cg =~ s/\n//g;
$ga = system("grep -o 'GA' sequence.txt | wc -l");
$ga =~ s/\n//g;
$at = system("grep -o 'AT' sequence.txt | wc -l");
$at =~ s/\n//g;
$tc = system("grep -o 'TC' sequence.txt | wc -l");
$tc =~ s/\n//g;
$a = system("grep -o 'A' sequence.txt | wc -l");
$a =~ s/\n//g;
$t = system("grep -o 'T' sequence.txt | wc -l");
$t =~ s/\n//g;
$g = system("grep -o 'G' sequence.txt | wc -l");
$g =~ s/\n//g;
$c = system("grep -o 'C' sequence.txt | wc -l");
$c =~ s/\n//g;
$pal = system("grep -o 'GCGC' sequence.txt | wc -l");#ABUNDANCIA HIP1
$pal =~ s/\n//g;
$EXP = (($gc*$cg*$ga*$at*$tc*$cg*$gc)/($c*$g*$a*$t*$c*$g))*($conteo-4+1);#ABUNDANCIA ESPERADA (EXP)

my $renglon = "$gc,$cg,$ga,$at,$tc,$a,$t,$g,$c,$conteo,$EXP,$pal";

open (MAR, ">>abundancia.csv") or die ("No puedo abrir abundancia.csv\n");
	print MAR ("sequence.txt\n\nGC,CG,GA,AT,TC,A,T,G,C,LEN,EXP,OBS\n");
	print MAR ("$renglon");
close (MAR);

close (LUZ);
