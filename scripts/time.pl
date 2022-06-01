#!/usr/local/bin/perl
 
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

my $lastdate = "$mday-$mon-$year\_$hour\:$min";
system("mkdir Archivos_$lastdate");


$genome_length = length(.fasta);
