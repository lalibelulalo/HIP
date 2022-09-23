#!/usr/local/bin/perl
use File::Copy;
# perl assembly_summary_sublist.pl assembly_summary.txt spp_list.txt outname
my $assembly_summary = $ARGV[0];
my $file_assembly = $assembly_summary;
$file_assembly =~ s/\.txt//g;
my $spp_list = $ARGV[1];
my $out_name = $ARGV[2];

print("$assembly_summary\n");
print("$spp_list\n");
print("$out_name\n");
system ("awk -F \"\\t\" '\$12==\"Complete Genome\" && \$11==\"latest\" {print \$0}' $assembly_summary >$file_assembly\_latest_complete.txt");

open (DIA, ">spp_list_2.tmp") or die ("No puedoabrir spp_list_2.tmp\n");
	open (KAR, "$spp_list") or die ("No puedoabrir $spp_list\n");
		while (my $spp = <KAR>){
			chomp ($spp);
			$spp =~ s/\)/\\\)/g;
			$spp =~ s/\(/\\\(/g;
			$spp =~ s/\'/\\047/g;
			print DIA ("$spp\n");
		}
	close (KAR);
close (DIA);

open (KAR, "spp_list_2.tmp") or die ("No puedoabrir spp_list_2.tmp\n");
	while (my $spp = <KAR>){
		chomp ($spp);
		system ("awk -F \"\\t\" '\$8==\"$spp\" {print \$0}' $file_assembly\_latest_complete.txt >> $file_assembly\_latest_complete_$out_name.txt");
		system ("awk -F \"\\t\" '\$8==\"$spp\" {print \$0}' $assembly_summary >> $file_assembly\_$out_name.txt");
	}
close (KAR);

unlink "spp_list_2.tmp";

