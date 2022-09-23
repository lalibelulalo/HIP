#!/usr/local/bin/perl

#my $spp = q(Acaryochloris sp. 'Moss Beach');
system ("awk -F \"\\t\" '\$8==\"Acaryochloris sp. \047Moss Beach\047\" {print \$0}' assembly_summary_latest_complete.txt");
#system ("awk -F \"\\t\" '\$8==\"Dolichospermum flos-aquae CCAP 1403/13F\" {print \$0}' assembly_summary_latest_complete.txt");
#system ("awk -F \"\\t\" '\$8==\"Acaryochloris marina MBIC11017\" {print \$0}' assembly_summary_latest_complete.txt");

#system ("awk -F \"\\t\" '\$12==\"Complete Genome\" && \$11==\"latest\" {print \$8\"_\"\$9\"_\"\$10}' $assembly_summary > species.txt");

#awk -F "\t" '$8=="Acaryochloris sp. \047Moss Beach\047" {print $0}' assembly_summary_latest_complete.txt
