#!/usr/local/bin/perl
use File::Copy;

my $assembly_summary = $ARGV[0]; ## Archivo assembly_summary.txt de refseq/genbank
my $out_name = $ARGV[1]; ## Nombre de salida para la carpeta

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(); ## Obtenemos la fecha para la carpeta de salida y tener un regisro
if ($min<10){
	$min = "0$min";
}
my $lastdate = "$mday-$mon-$year\_$hour\:$min";
system("mkdir NCBI_$out_name\_Genome_Files_$lastdate"); ## Creamos la carpeta con la fecha y el nombre de salida proporcionado

system("awk -F \"\\t\" '{print \$8\"_\"\$9\"_\"\$10}' $assembly_summary > species.txt"); ## Extraemos las especies y las graduamos en el archivo species.txt
system("awk -F \"\\t\" '{print \$20}' $assembly_summary > ftp_dirs.txt"); ## Extraemos las direcciones ftp y las graduamos en el archivo ftp_dirs.txt
system("sed -i 's/\\s/_/g' species.txt"); ## quitamos espacios y caracteres especiales para poder procesar las especies
system("sed -i 's/\\./_/g' species.txt");
system("sed -i 's/:/_/g' species.txt");
system("sed -i 's/strain=//g' species.txt");
system("sed -i 's/\\s/_/g' species.txt");
system("sed -i 's/__/_/g' species.txt");
system("sed -i 's/\\//-/g' species.txt");
system("sed -i 's/(/_/g' species.txt");
system("sed -i 's/)//g' species.txt");
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
	print $SppFtpFile join(",", $columns1[1-1], $columns2[1-1]),"\n";
}

close($SppFile);
close($FtpFile);
close($SppFtpFile);
=pod
assembly_structure/                         2022-07-30 10:03    -   
assembly_report.txt                         2020-02-04 15:34   11K  
assembly_stats.txt                          2022-01-06 09:40  4.8K  
cds_from_genomic.fna.gz                     2022-07-30 10:03  1.3M  
feature_count.txt.gz                        2022-07-30 10:03  263   
feature_table.txt.gz                        2022-07-30 10:03  177K  
genomic.fna.gz                              2018-09-04 05:12  1.2M  
genomic.gbff.gz                             2022-07-30 10:03  3.0M  
genomic.gff.gz                              2022-07-30 10:03  310K  
genomic.gtf.gz                              2022-07-30 10:03  366K  
genomic_gaps.txt.gz                         2020-02-04 15:34  484   
protein.faa.gz                              2022-07-30 10:03  809K  
protein.gpff.gz                             2022-07-30 10:03  2.0M  
rna_from_genomic.fna.gz                     2021-01-06 01:33  4.5K  
translated_cds.faa.gz                       2022-07-30 10:03  926K
=cut
open (KAR, "genome_dirs.txt") or die ("No puedoabrir genome_dirs.txt\n");
	while (my $line = <KAR>){
		chomp ($line);
		if ($line =~ /(.+),(.+)/){
			$spp=$1;
			$link=$2;
			}
		if ($link =~ /(\w\w\w_.+)/){
				$acc=$1;
				}
		mkdir $spp;
		#print ("$spp\t$acc\t$link\n");
		system ("wget $link/$acc\_protein.faa.gz");
		move ("$acc\_protein.faa.gz", "./$spp");
		
		system ("wget $link/$acc\_genomic.fna.gz ");
		move ("$acc\_genomic.fna.gz", "./$spp");
		
		system ("wget $link/$acc\_genomic.gbff.gz");
		move ("$acc\_genomic.gbff.gz", "./$spp");
		
		system ("wget $link/$acc\_genomic.gff.gz");
		move ("$acc\_genomic.gff.gz", "./$spp");
		
		system ("wget $link/$acc\_genomic.gtf.gz");
		move ("$acc\_genomic.gtf.gz", "./$spp");
		
		system ("wget $link/$acc\_cds_from_genomic.fna.gz");
		move ("$acc\_cds_from_genomic.fna.gz", "./$spp");
		#----------------------------------
		system ("gzip -d ./$spp/*.gz");
		#----------------------------------	
		rename ("./$spp/$acc\_protein.faa", "./$spp/$spp\_protein.faa");
		rename ("./$spp/$acc\_genomic.fna", "./$spp/$spp\_genomic.fna");
		rename ("./$spp/$acc\_genomic.gbff", "./$spp/$spp\_genomic.gbff");
		rename ("./$spp/$acc\_genomic.gff", "./$spp/$spp\_genomic.gff");
		rename ("./$spp/$acc\_genomic.gtf", "./$spp/$spp\_genomic.gtf");
		rename ("./$spp/$acc\_cds_from_genomic.fna", "./$spp/$spp\_cds_from_genomic.fna");
			
		#----------------------------------
		move ("./$spp", "./NCBI_$out_name\_Genome_Files_$lastdate/$spp");
		
			
		}
close (KAR);

unlink "species.txt";
unlink "ftp_dirs.txt";
unlink "genome_dirs.txt";
