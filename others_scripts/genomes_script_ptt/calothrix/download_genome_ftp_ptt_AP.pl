#!/usr/local/bin/perl

		## 1 - NO REGEX
		## Find: '
		## Replace: \'

		## 2 - REGEX
		## In GBKtoPTT.pl:
		## Find: (\t *+)(.+)
		## Replace: print MAR ('\1\2');

		## 3 - REGEX
		## Find: \t
		## Replace: nothing

		## 4 - NO REGEX
		## Find: ');\n
		## Replace: ',"\\n");\n

open (MAR, ">GBKtoPTT.pl") or die ("GBKtoPTT.pl\n");
	print MAR ('#!/usr/bin/env perl',"\n");
	print MAR ('use strict;',"\n");
	print MAR ('use Bio::SeqIO;',"\n");
	print MAR ("\n");
	print MAR ('# This script takes a GenBank file as input, and produces a',"\n");
	print MAR ('# NCBI PTT file (protein table) as output. A PTT file is',"\n");
	print MAR ('# a line based, tab separated format with fixed column types.',"\n");
	print MAR ('#',"\n");
	print MAR ('my $gbk = Bio::SeqIO->new(-fh=>\*STDIN, -format=>\'genbank\');',"\n");
	print MAR ('my $seq = $gbk->next_seq;',"\n");
	print MAR ('my @cds = grep { $_->primary_tag eq \'CDS\' } $seq->get_SeqFeatures;',"\n");
	print MAR ("\n");
	print MAR ('print $seq->description, " - 0..",$seq->length,"\r\n";',"\n");
	print MAR ('print scalar(@cds)," proteins\r\n";',"\n");
	print MAR ('print join("\t", qw(Location Strand Length PID Gene Synonym Code COG ',"\n");
	print MAR ('Product)),"\r\n";',"\n");
	print MAR ("\n");
	print MAR ('for my $f (@cds) {',"\n");
	print MAR ('   my $gi = \'-\';',"\n");
	print MAR ('   $gi = $1 if tag($f, \'db_xref\') =~ m/\bGI:(\d+)\b/;',"\n");
	print MAR ('   my $cog = \'-\';',"\n");
	print MAR ('   $cog = $1 if tag($f, \'product\') =~ m/^(COG\S+)/;',"\n");
	print MAR ('   my @col = (',"\n");
	print MAR ('     $f->start.\'..\'.$f->end,',"\n");
	print MAR ('     $f->strand >= 0 ? \'+\' : \'-\',',"\n");
	print MAR ('     ($f->length/3)-1,',"\n");
	print MAR ('     $gi,',"\n");
	print MAR ('     tag($f, \'gene\'),',"\n");
	print MAR ('     tag($f, \'locus_tag\'),',"\n");
	print MAR ('     $cog,',"\n");
	print MAR ('     tag($f, \'product\'),',"\n");
	print MAR ('   );',"\n");
	print MAR ('   print join("\t", @col), "\r\n";',"\n");
	print MAR ('}',"\n");
	print MAR ("\n");
	print MAR ('sub tag {',"\n");
	print MAR ('   my($f, $tag) = @_;',"\n");
	print MAR ('   return \'-\' unless $f->has_tag($tag);',"\n");
	print MAR ('   return join(\' \', $f->get_tag_values($tag));',"\n");
	print MAR ('}',"\n");
close (MAR);

mkdir("PGAP_files");
mkdir("genbank");
system ("mv genbank PGAP_files");
mkdir("PGAP-x_files");
open (LAL, "genomes_list.txt") or die ("No puedoabrir genomes_list\n");
	while (my $genome = <LAL>){
		$genome =~ s/\n//g;
		
		my $A;
		my $B;
		my $C;
		my $D;
		my $E;
		my $feature_tab;
		my $protein_faa;
		my $genomic_fna;
		my $cds_genomic;
		my $transla_cds;
		my $genbank;
		my $especie;

		if ($genome =~ /([A-Z0-9]+)_([A-Z0-9][A-Z0-9][A-Z0-9])([A-Z0-9][A-Z0-9][A-Z0-9])([A-Z0-9][A-Z0-9][A-Z0-9])\.(\d+)_(\w+)\t(.+)/){
			$A=$1;	# AAA
			$B=$2;	# 111
			$C=$3;	# 222
			$D=$4;	# 333
			$E=$5;	# 4
			$F=$6;	# XXXXXX
			$G=$7;	# ESPECIE => AAA_111222333.4_XXXXXX	ESPECIE
			$especie="$G";
			$genome="$A\_$B$C$D.$E\_$F";
			print("GENOMA: $genome\n");
			print("ESPECIE: $especie\n");

			$feature_tab="https://ftp.ncbi.nlm.nih.gov/genomes/all/$1/$2/$3/$4/$genome/$genome\_feature_table.txt.gz";
			$protein_faa="https://ftp.ncbi.nlm.nih.gov/genomes/all/$1/$2/$3/$4/$genome/$genome\_protein.faa.gz";
			$genomic_fna="https://ftp.ncbi.nlm.nih.gov/genomes/all/$1/$2/$3/$4/$genome/$genome\_genomic.fna.gz";
			$cds_genomic="https://ftp.ncbi.nlm.nih.gov/genomes/all/$1/$2/$3/$4/$genome/$genome\_cds_from_genomic.fna.gz";
			$transla_cds="https://ftp.ncbi.nlm.nih.gov/genomes/all/$1/$2/$3/$4/$genome/$genome\_translated_cds.faa.gz";
			$genbank_gbf="https://ftp.ncbi.nlm.nih.gov/genomes/all/$1/$2/$3/$4/$genome/$genome\_genomic.gbff.gz";

			system ("curl $feature_tab >$especie\_$genome\_feature_table.txt.gz");
			system ("curl $protein_faa >$especie\_$genome\_protein.faa.gz");
			system ("curl $cds_genomic >$especie\_$genome\_cds_from_genomic.fna.gz");
			system ("curl $genbank_gbf >$especie\_$genome\_genomic.gbff.gz");
			system ("curl $genbank_gbf >$especie\_$genome\_genomic.fna.gz");
			system ("gzip -d $especie\_$genome\_genomic.gbff.gz");
			system ("gzip -d $especie\_$genome\_genomic.fna.gz");
			system ("perl GBKtoPTT.pl < $especie\_$genome\_genomic.gbff > $especie\_$genome\_genomic.ptt");
			mkdir($especie);
			
			system ("cp *_feature_table.txt.gz PGAP_files/
				cp *_protein.faa.gz PGAP_files/
				cp *_cds_from_genomic.fna.gz PGAP_files/
				cp $especie\_$genome\_genomic.gbff PGAP_files/genbank/
				cp $especie\_$genome\_genomic.fna PGAP-x_files/
				cp *_genomic.ptt PGAP-x_files");
				
			system ("mv PGAP_files/genbank/$especie\_$genome\_genomic.gbff PGAP_files/genbank/$especie\_$genome\_genomic.gb");
			system ("mv *.gz $especie/
				mv *.gbff $especie/
				mv *.ptt $especie/
				mv *.fna $especie");

				}
		}
close (LAL);
system ("rm GBKtoPTT.pl");
			
			
			
			
			
			
			
			
=pod
$genome\_assembly_report.txt
$genome\_ASM73489v2_assembly_stats.txt 
$genome\_ASM73489v2_cds_from_genomic.fna.gz
$genome\_ASM73489v2_feature_count.txt.gz
$genome\_ASM73489v2_feature_table.txt.gz
$genome\_ASM73489v2_genomic.fna.gz
$genome\_ASM73489v2_genomic.gbff.gz
$genome\_ASM73489v2_genomic.gff.gz
$genome\_ASM73489v2_genomic.gtf.gz
$genome\_ASM73489v2_protein.faa.gz
$genome\_ASM73489v2_protein.gpff.gz
$genome\_ASM73489v2_rna_from_genomic.fna.gz
$genome\_ASM73489v2_translated_cds.faa.gz
README.txt
annotation_hashes.txt
assembly_status.txt 
md5checksums.txt
=cut
