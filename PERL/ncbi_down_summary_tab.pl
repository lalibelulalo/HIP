#!/usr/local/bin/perl
use File::Copy;

my $assembly_summary = $ARGV[0]; ## Archivo assembly_summary.txt de refseq/genbank
my $acc_spp_table = $ARGV[1]; ## tabla de ncbi 
my $out_name = $ARGV[2]; ## Nombre de salida para la carpeta

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(); ## Obtenemos la fecha para la carpeta de salida y tener un regisro
if ($min<10){
	$min = "0$min";
}
my $lastdate = "$mday-$mon-$year\_$hour\:$min";
system("mkdir NCBI_$out_name\_Genome_Files_$lastdate"); ## Creamos la carpeta con la fecha y el nombre de salida proporcionado

system("awk -F \"\\t\" '{print \$8\"_\"\$9\"_\"\$10}' $assembly_summary > species.txt"); ## Extraemos las especies y las guardamos en el archivo species.txt
system("awk -F \"\\t\" '{print \$1\"\\t\"\$20}' $assembly_summary > ftp_dirs.txt"); ## Extraemos las direcciones ftp y las graduamos en el archivo ftp_dirs.txt
system("sed -i 's/\\s/_/g' species.txt"); ## quitamos espacios y caracteres especiales para poder procesar las especies
system("sed -i 's/\\./_/g' species.txt");
system("sed -i 's/:/_/g' species.txt");
system("sed -i 's/;/_/g' species.txt");
system("sed -i 's/strain=//g' species.txt");
system("sed -i 's/\\s/_/g' species.txt");
system("sed -i 's/__/_/g' species.txt");
system("sed -i 's/\\//-/g' species.txt");
system("sed -i 's/(/_/g' species.txt");
system("sed -i 's/\\[/_/g' species.txt");
system("sed -i 's/)//g' species.txt");
system("sed -i 's/=//g' species.txt");
system("sed -i 's/\\]/_/g' species.txt");
system("sed -i 's/\\x27/_/g' species.txt");
system("sed -i 's/#//g' species.txt");
system("sed -i 's/__/_/g' species.txt");

open my $SppFile, '<', "species.txt" or die qq{Failed to open "species.txt" for writing: $!};
open my $FtpFile, '<', "ftp_dirs.txt" or die qq{Failed to open "ftp_dirs.txt" for writing: $!};
open my $SppFtpFile, '>', "genome_dirs.txt" or die qq{Failed to open "genome_dirs.txt" for writing: $!};

# Uno especies y links
while(my $l1 = <$SppFile>){ 
	my $l2 = <$FtpFile>;
	chomp $l1;
	chomp $l2;
	my @columns1 = split(/ /, $l1);
	my @columns2 = split(/ /, $l2);
	print $SppFtpFile join("\t", $columns1[1-1], $columns2[1-1]),"\n";
}

close($SppFile);
close($FtpFile);
close($SppFtpFile);

my $contador = 0;
open (LUZ, "genome_dirs.txt") or die ("No puedo abrir genome_dirs.txt");
	while (my $line1 = <LUZ>){
		chomp ($line1);
		my @assembly_sum = split (/\t/, $line1);
		$acc_sum = @assembly_sum[1];
		$spp = @assembly_sum[0];
		$link = @assembly_sum[2];
		if ($link =~ /(\w\w\w_.+)/){
				$acc_sum2=$1;
				}
		
		open (DIA, "$acc_spp_table") or die ("No puedo abrir $acc_spp_table\n");
			while (my $line2 = <DIA>){
				chomp ($line2);
				my @tab_sum = split (/\t/, $line2);
				$acc_tab = @tab_sum[0];
				
				if ($acc_sum eq $acc_tab){
					mkdir $spp;
					$contador++;
					print ("$spp\t$acc_sum2\t$link\n");
					system ("wget $link/$acc_sum2\_protein.faa.gz");
					move ("$acc_sum2\_protein.faa.gz", "./$spp");
					
					system ("wget $link/$acc_sum2\_genomic.fna.gz ");
					move ("$acc_sum2\_genomic.fna.gz", "./$spp");
					
					system ("wget $link/$acc_sum2\_genomic.gbff.gz");
					move ("$acc_sum2\_genomic.gbff.gz", "./$spp");
					
					system ("wget $link/$acc_sum2\_genomic.gff.gz");
					move ("$acc_sum2\_genomic.gff.gz", "./$spp");
					
					system ("wget $link/$acc_sum2\_genomic.gtf.gz");
					move ("$acc_sum2\_genomic.gtf.gz", "./$spp");
					
					system ("wget $link/$acc_sum2\_cds_from_genomic.fna.gz");
					move ("$acc_sum2\_cds_from_genomic.fna.gz", "./$spp");
					#----------------------------------
					system ("gzip -d ./$spp/*.gz");
					#----------------------------------	
					rename ("./$spp/$acc_sum2\_protein.faa", "./$spp/$spp\_protein.faa");
					rename ("./$spp/$acc_sum2\_genomic.fna", "./$spp/$spp\_genomic.fna");
					rename ("./$spp/$acc_sum2\_genomic.gbff", "./$spp/$spp\_genomic.gbff");
					rename ("./$spp/$acc_sum2\_genomic.gff", "./$spp/$spp\_genomic.gff");
					rename ("./$spp/$acc_sum2\_genomic.gtf", "./$spp/$spp\_genomic.gtf");
					rename ("./$spp/$acc_sum2\_cds_from_genomic.fna", "./$spp/$spp\_cds_from_genomic.fna");
						
					#----------------------------------
					move ("./$spp", "./NCBI_$out_name\_Genome_Files_$lastdate/$spp");

				}

				
			}
		close (DIA);
					}
close (LUZ);


unlink "species.txt";
unlink "ftp_dirs.txt";
unlink "genome_dirs.txt";
