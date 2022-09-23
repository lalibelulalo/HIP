#!/usr/local/bin/perl
use File::Copy;

my $ftp_dirs = $ARGV[0];
my $spp_list = $ARGV[1];
my $out_name = $ARGV[2];

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
if ($min<10){
	$min = "0$min";
}
my $lastdate = "$mday-$mon-$year\_$hour\:$min";
system("mkdir NCBI_$out_name\_Genome_Files_$lastdate");

open (DIA, ">spp_list.tmp") or die ("No puedoabrir spp_list.tmp\n");
	open (LAU, "$spp_list") or die ("No puedoabrir $spp_list\n");
		while (my $line = <LAU>){
			chomp ($line);
			print DIA ("$line\n");
			}
	close (LAU);
close (DIA);

system("sed -i 's/\\s/_/g' spp_list.tmp");
system("sed -i 's/\\./_/g' spp_list.tmp");
system("sed -i 's/:/_/g' spp_list.tmp");
system("sed -i 's/strain=//g' spp_list.tmp");
system("sed -i 's/\\s/_/g' spp_list.tmp");
system("sed -i 's/__/_/g' spp_list.tmp");
system("sed -i 's/\\//-/g' spp_list.tmp");
system("sed -i 's/(/_/g' spp_list.tmp");
system("sed -i 's/)//g' spp_list.tmp");
system("sed -i 's/#//g' spp_list.tmp");
system("sed -i 's/__/_/g' spp_list.tmp");

open my $input1, '<', "spp_list.tmp" or die qq{Failed to open "spp_list.tmp" for writing: $!};
open my $input2, '<', "$ftp_dirs" or die qq{Failed to open "$ftp_dirs" for writing: $!};
open my $outfile, '>', "genome_dirs.txt" or die qq{Failed to open "genome_dirs.txt" for writing: $!};

# Uno especies y links
while(my $l1 = <$input1>){ 
	my $l2 = <$input2>;
	chomp $l1;
	chomp $l2;
	my @columns1 = split(/ /, $l1);
	my @columns2 = split(/ /, $l2);
	print $outfile join(",", $columns1[1-1], $columns2[1-1]),"\n";
}

close($input1);
close($input2);
close($outfile);
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

unlink "$ftp_dirs";
unlink "genome_dirs.txt";
unlink "spp_list.tmp";
