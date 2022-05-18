Introduction
#############

¿Qué es PhaME?
==============

PhaME o Phylogenetic and Molecular Evolution (PhaME) es una herramienta de análisis que realiza análisis filogenéticos y de evolución molecular.

Dada una referencia, PhaME extrae SNP de genomas completos, borradores de genomas y/o lecturas, utiliza alineación de secuencias múltiples de SNP para construir un árbol filogenético y proporciona análisis evolutivos (genes bajo selección positiva) utilizando CDS SNP.


PhaME: Under the hood
======================

PhaME’s input consists of a set of genomes in fasta and/or fastq formats, and corresponding annotation files in gff3 format if downstream analyses will include coding regions. A detailed step by step explanation of how PhaME analyzes genomes is provided below.

Selecting Reference genome:
-----------------------------
Since PhaME is a reference genome-based tool where all input genomes and metagenomes are aligned against a reference, the first step of PhaME’s analysis is selecting a reference genome. Given a set of genomes (in a folder listed under the “refdir” parameter of the control file), the reference genome can be selected using one of three options: option 1- a random genome is selected from the provided set of genomes; option 2- a specific genome is selected from the set via input from the user; option 3- the MinHash distance is calculated between all genomes provided (complete genomes, draft genomes, and raw reads) to determine which reference genome has the shortest average distance to all of the other genomes. MinHash distances are calculated using its implementation in BBMap [Bushnell]_.

Self-nucmerization to remove repeats from reference genomes:
---------------------------------------------------------------
The genome alignment portion of PhaME is built on the tool nucmer2 for alignments of genomes in FASTA format. Each genome included is first aligned with itself using nucmer, called self-nucmerization, and then aligned regions called repeats are removed from the genomes for downstream analyses. The following nucmer command is used for the self-nucmerization step: 

::

    nucmer --maxmatch --nosimplify --prefix=seq_seq ref_genomeA.fasta ref_genomeA.fasta 

::

The option --maxmatch, which reports all matches, is used to ensure that all possible alignments are reported for maximal removal of repeats. 

Genome Alignments
--------------------------------
All genomes that are in `refdir` are first aligned against the reference genome (see section 1) that have had its repeats removed (section 2). Likewise, incomplete genomes or contigs, the ones that are listed in the `workdir` with extension `.contig` are also aligned against the reference genome using `nucmer`. For aligning genomes in fasta format against each other, following command, same as the previous step for nucmer alignment is used:

::

    nucmer --maxmatch refgenome.fasta genome.fasta

::

All other options in nucmer alignments are kept at default, some of the important ones are listed below:

::

   -b|breaklen     Set the distance an alignment extension will attempt to
                    extend poor scoring regions before giving up (default 200)
    -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
    -D|diagdiff     Set the maximum diagonal difference between two adjacent
                    anchors in a cluster (default 5)
    -d|diagfactor   Set the maximum diagonal difference between two adjacent
                    anchors in a cluster as a differential fraction of the gap
                    length (default 0.12)
    --[no]extend    Toggle the cluster extension step (default --extend)
    -g|maxgap       Set the maximum gap between two adjacent matches in a
                    cluster (default 90)
    -l|minmatch     Set the minimum length of a single match (default 20)

::

Also, `nucmer` only aligns `"A"` `"T"` `"C"` `"G"`, all other characters are ignored. So, if there are `"N"`s in the provided genomes, thse positions are not included in the alignment.

*Note*: If an analysis requires running multiple iterations of PhaME on a same set of dataset or a subset of dataset, one doesn't need to perform the alignment over and over again. PhaME provides an option where it can keep all possible pairwise alignment of genomes from `refdir` for future analyses. All the steps mentioned in this section are the same, except that all vs. all alignment is performed compared to just one reference. By doing all vs. all alignment one can also test the effect on their analyses with different reference genomes.

Mapping of raw reads to reference genome
-------------------------------------------
Currently, PhaME only processes short, raw reads from Illumina. If raw reads, single or paired end, are included in the analyses, they are mapped to the reference genome using either bowtie2 or BWA based on users’ input. For reads mapping to the reference genome, the following commands are used:

First, it builds database from the reference genome.
::

    bowtie2-build refgenome refgenome

::
or, if BWA was chosen as the preferred aligner:

::

    bwa index refgenome

::

The raw reads are then mapped to the reference genomne using one of the following commands:

For bowtie2 and paired reads:

::

    bowtie2 -a -x $refgenome -1 read1 -2 read2 -S paired.sam`;

::
The option `-a` reports all possible alignments.

For bowtie2 and single end reads:

::

    bowtie2 -a -x $refgenome -U read -S single.sam`;

::

For BWA and paired reads:

::

    bwa mem refgenome read1 read2 | samtools view -ubS -| samtools sort -T tmp_folder -O BAM -o paired.bam

::

For BWA and single end reads:

::

    bwa mem refgenome read |samtools view -ubS - | samtools sort -T tmp_folder -O BAM -o single.bam

::


Filtering genome alignments
------------------------------
Genome alignment produced using `nucmer` are filtered using `delta-filter` to only keep 1 to 1 alignments allowing for rearrangements. This filtering step is produced for all `nucmer` alignments.

::

    delta-filter -1 genome.delta > genome.snpfilter

::


Calling SNPs from genome alignments
--------------------------------------
The pairwise `nucmer` alignments are then parsed to produce a SNP table using `show-snps`.

::

    show-snps -CT genome.snpfilter > genome.snps

::

Here, option C and T specifies not to report SNPs from ambiguous alignments and report the output in tab delimited file respectively.

Reporting nucmer alignments
------------------------------

Each alignments are further parse to produce a tab delimited file that has information on regions and %ID of their alignments.
::

    show-coords -clTr genome.snpfilter > genome.coords

::

The parameter flag -clTr implies different headers to be reported in the report.

::

-c          Include percent coverage information in the output
-l          Include the sequence length information in the output
-r          Sort output lines by reference IDs and coordinates
-T          Switch output to tab-delimited format

::

Calling SNPs from read mapping
---------------------------------
`bcftools mpileup` is used for calling SNPs from read mapping results (bam file) of every genomes represented by raw reads. Maximum depth is set to 1000000 for both SNP and indel calling and minimum gaps for calling an indel is set to 3. The output vcf file is then used to call SNPs using `bcftools call` where ploidy is specified as `1` if its a haploid or bacterial genome, else it is called using default parameter. Furthermore, based on the user specified parameter in the control file, SNPs are further filtered based on percentage of SNPs. Here are the snippets of commmand that are run as part of this. All of them result in a vcf file.

::

    bcftools mpileup -d 1000000 -L 1000000 -m 3 -Ov -f $refgenome $bam_output | bcftools call --ploidy 1 -cO b > $bcf_output;
    bcftools view -v snps,indels,mnps,ref,bnd,other -Ov $bcf_output | vcfutils.pl varFilter -a$min_alt_bases -d$min_depth -D$max_depth > $vcf_output`;
    bcftools filter -i '(DP4[0]+DP4[1])==0 || (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) > $snp_filter' $vcf_output > $vcf_filtered`

::


Calculating core genome alignments
----------------------------------
As a first step in calculating the core genome, all alignments to reference are checked for linear coverage to assure the proportion of reference genome that was used in the alignment. If its lower than the threshold cutoff (default = 0.6) set in control file, that genome will be removed from further analyses. Then rest of the pairwise alignments that are either in vcf format or nucmer formats are then collated to calculate a core genome. Only the alignment position that are 100% conserved are kept, all other positions are removed from the final core genome alignment. PhaME produces multiple alignment files corresponding to core genome such as the one that has only the variant sites (`_all_snp_alignment.fna`), has variant and invariant sites (`all_alignment.fna`), and the ones that have SNPs from only the coding region (`_cds_snp_alignment.fna`). The coding region SNP alignment requires a GFF formatted annotation file.


Reconstructing core genome phylogeny
-------------------------------------
PhaME provides multiple tools (RAxML [Stamatakis 2014]_, FastTree [Price 2010]_, and IQ-Tree [Nguyen 2015]_) to reconstruct phylogeny from one core genome alignments that have invariant sites. If RAxML or FastTree option is chosen, users cannot modify the models as they are pre-selected. RAxML trees are reconstructed using GTRGAMMAI models that "GTR + Optimization of substitution rates + GAMMA model of rate heterogeneity (alpha parameter will be estimated)" with `I` but with estimate for invariable sites. FastTree uses GTR model only. IQ-TREE is run using option `-m TEST` that searches for the best model that fits the data before reconstructing the phylogeny. RAxML is the only option that is currently available that can also calculate the bootstraps.

Selecting genes for molecular evolutionary analyses
-------------------------------------------------------
To perform selection analyses using PAML or HyPhy, codon alignments of genes are required. Based on the position of SNPs in the reference genome, if a SNP is within a coding region and if that coding region does not have a gap, they are extracted from the core genome alignment. The nucleotide sequences of the genes are translated to protein sequences, aligned using the program mafft 8, and then reverse translated back to nucleotide using the Perl code pal2nal.pl from http://www.bork.embl.de/pal2nal/.

Molecular Evoluationary analyses
------------------------------------

The set of gene alignments are used for molecular evolutionary analyses using either PAML [Yang 2007]_ or HyPhy. Both packages can test for the presence of positively selected sites and lineages by allowing the dN/dS ratio (ω) to vary among sites and lineages. The adaptive branch-site REL test for episodic diversification (aBSREL) model in the HyPhy package is used to detect instances of episodic diversifying and positive selection. If PAML is selected, the M1a-M2a and M7- M8 nested models are implemented. In the latter case, the likelihood ratio test between the null models (M1a and M8) and the alternative model (M2a and M7) at a significance cutoff of 5% provides information on how the genes are evolving. The results for each gene are then summarized in a table containing information on whether the gene is evolving under positive, neutral, or negative selection, along with p-values. HyPhy is run with a model, which specifically looks for sign of positive selection in given sets of genes. The analysis produces a list of JSON files corresponding to each gene which can be uploaded to vision.hyphy.org/absrel for further analysis. We opted to provide PAML as an option, however we recommend using HyPhy for large projects due to its speed and concise output. 


References
--------------
.. [Yang 2007] Yang Z: PAML 4: phylogenetic analysis by maximum likelihood. Mol Biol Evol 2007, 24:1586-1591.
.. [Pond 2005] Pond SL, Frost SD, Muse SV: HyPhy: hypothesis testing using phylogenies. Bioinformatics 2005, 21:676-679.
.. [Kurtz 2004] Kurtz S, Phillippy A, Delcher AL, Smoot M, Shumway M, Antonescu C, Salzberg SL: Versatile and open software for comparing large genomes. Genome Biol 2004, 5:R12.
.. [Bushnell] Bushnell B: BBMap. 37.66 edition. sourceforge.net/projects/bbmap/.
.. [Stamatakis 2014] Stamatakis A: RAxML version 8: a tool for phylogenetic analysis and post- analysis of large phylogenies. Bioinformatics 2014, 30:1312-1313.
.. [Price 2010] Price MN, Dehal PS, Arkin AP: FastTree 2--approximately maximum- likelihood trees for large alignments. PLoS One 2010, 5:e9490.
.. [Nguyen 2015] Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ: IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Mol Biol Evol 2015, 32:268-274.
