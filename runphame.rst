How to run PhaME?
#################
After installation, PhaME only requires a control file along with your input files.

.. code-block:: console

    $ phame phame.ctl


What is in control file?
========================
.. _control_file:
    .. In a control file, parameters and input folders are specified.

    .. Here is how a control file looks like with the description of their options.


::

       refdir = test/data/ebola_ref  # directory where reference (Complete) files are located
      workdir = test/workdirs/ebola_complete # directory where contigs/reads files are located and output is stored

    reference = 2  # 0:pick a random reference from refdir; 1:use given reference; 2: use ANI based reference
      reffile = KJ660347.fasta  # reference filename when option 1 is chosen

      project = t4  # main alignment file name

      cdsSNPS = 1  # 0:no cds SNPS; 1:cds SNPs, divides SNPs into coding and non-coding sequences, gff file is required

      buildSNPdb = 0 # 0: only align to reference 1: build SNP database of all complete genomes from refdir

    FirstTime = 1  # 1:yes; 2:update existing SNP alignment, only works when buildSNPdb is used first time to build DB

         data = 0  # *See below 0:only complete(F); 1:only contig(C); 2:only reads(R); 
                   # 3:combination F+C; 4:combination F+R; 5:combination C+R; 
                   # 6:combination F+C+R; 7:realignment  *See below 
        reads = 2  # 1: single reads; 2: paired reads; 3: both types present;

         tree = 1  # 0:no tree; 1:use FastTree; 2:use RAxML; 3:use both;
    bootstrap = 0  # 0:no; 1:yes;  # Run bootstrapping  *See below
            N = 100  # Number of bootstraps to run *See below    
  
    PosSelect = 1  # 0:No; 1:use PAML; 2:use HyPhy; 3:use both # these analysis need gff file to parse genomes to genes

         code = 0  # 0:Bacteria; 1:Virus; 2: Eukarya # Bacteria and Virus sets ploidy to haploid

        clean = 0  # 0:no clean; 1:clean # remove intermediate and temp files after analysis is complete

      threads = 2  # Number of threads to use

       cutoff = 0.1  # Linear alignment (LA) coverage against reference - ignores SNPs from organism that have lower cutoff.

       * When using data option 1,2,5 need a complete reference to align/map to.
       * Use data option 7 when need to extract SNPs using a sublist of already aligned genomes.


It is a simple text file similar to the ones used PAML analysis.

What do parameters in control file stand for?
==============================================

1. *refdir*
A directory with reference genomes (complete genomes) and their annotation file in gff format (optional). Each file should represent a genome and have following extensions. The path to the folder should be either absolute or relative to the location of the control file. Please avoid filenames that have multiple `.` or has special characters like `:` in their name.
    - `*`.fasta
    - `*`.fna
    - `*`.fa
    - `*`.gff  (optional: to analyze Coding region SNPs of a selected reference file)    

For example, a typical reference folder with reference genomes look like this:

::

    $ ls ref/
    GCA_000006925_2_ASM692v2_genomic.fna   GCA_000017745_1_ASM1774v1_genomic.fna  GCA_000026245_1_ASM2624v1_genomic.fna   GCA_000227625_1_ASM22762v1_genomic.fna
    GCA_000007405_1_ASM740v1_genomic.fna   GCA_000017765_1_ASM1776v1_genomic.fna  GCA_000026265_1_ASM2626v1_genomic.fna   GCA_000245515_1_ASM24551v1_genomic.fna
    GCA_000008865_1_ASM886v1_genomic.fna   GCA_000017985_1_ASM1798v1_genomic.fna  GCA_000026265_1_ASM2626v1_genomic.gff   GCA_000257275_1_ASM25727v1_genomic.fna


Each of these files represent one genome. Each genome may have multiple sequences representing multiple replicons or contigs, but are all part of one genome. `gff` files corresponding to a genome must have the same exact name and in the same folder, just different extension. For example, `gff` file for genome `GCA_000006925_2_ASM692v2_genomic.fna` is `GCA_000006925_2_ASM692v2_genomic.gff`.

2. *workdir*
This is the folder where intermediate and final results of analysis are stored. The path to the folder should be either absolute or relative to the location of the control file. Additionally, if the analysis includes incomplete genomes or contig files and raw reads, they must be in this folder. Contigs file must have following extensions to be recognised as contig file.
     - `*`.contig
     - `*`.contigs

    For example, a working directory with contigs folder look like this:

::

    $ ls workdir/*.contig\
    workdir/GCA_000155105_1_ASM15510v1_genomic.contig  workdir/GCA_000968895_2_ASM96889v2_genomic.contig   workdir/GCA_001514825_1_ASM151482v1_genomic.contig
    workdir/GCA_000190495_1_ASM19049v1_genomic.contig  workdir/GCA_000968905.2_ASM96890v2_genomic.contig   workdir/GCA_001514845_1_ASM151484v1_genomic.contig
    workdir/GCA_000191665_1_ecmda7_genomic.contig      workdir/GCA_001471755_1_ASM147175v1_genomic.contig  workdir/GCA_001514865_1_ASM151486v1_genomic.contig


If the analysis includes reads, they must be in `workdir` as well and decompressed. If reads are paired, they must have same file name at the beginning of the name and `R1` and `R2` at the end of the name and needs to have `.fastq` as their extension (`*_`R1.fastq `*_`R2.fastq). Any file that have `*.fastq` as their extension but dont have paired reads will be treated as single reads. For example, a working folder with paired raw read files loole like this:

::

    $ ls *.fastq
    GGB_SRR2000383_QC_trimmed_R1.fastq  GGB_SRR2000383_QC_trimmed_R2.fastq  GGC_SRR2164314_QC_trimmed_R1.fastq  GGC_SRR2164314_QC_trimmed_R2.fastq


3. *reference*
    This is where you specify how do you want to pick your reference genome. The available options are:
        - 0: randomly pick a genome from `refdir` folder as the reference genome.
        - 1: use the specified genome as the reference. Genome's filename is specified in the `reffile` option.
        - 2: picks a `mid point` genome based on the  Average Nucleotide Identity (ANI) among all genomes. It uses mash (implemented in BBMap) to calculate ANI.

4. *reffile*
    This is where you specify the reference genome, if option 1 is picked in previous option. File name of the genome is written here and the program will look for that file in `reffile` folder. For example, `KJ660347.fasta` in the control file example above is found in the `reffile` folder.

5. *project*
    The name of the project. All the important downstream output filenames will have the specified project name as their prefix.

6. *cdsSNPS*
    This option allows users to parse SNPs based on their position into coding and non-coding sequences. It can be turned ON (0) or OFF (1). If turned ON, the picked reference genome must have a corresponding gff file. This option is automatically turned ON, if Molecular evolutionary analyses is turned ON (see below).

7. *buildSNPdb*
    This option will turn ON (1) or OFF (0) database creation, which is essentially all possible pairwise alignment of all genomes in `refdir`. Turning this ON will significantly increase the runtime.

8. *FirstTime*
    This options default is 1, which reruns everything. The option 2, which only recalculates the SNP matrix only works when SNP database is turned ON in previous step.

9. *data*
    Select the appropriate option based on the type of data that was included in the analysis. 
        - 0: only full/complete(F);

            + Select this option if you only have full/complete genomes or you only want to analyze these genomes from the dataset. Full/COmplete genomes are the ones that are found in *refdir*.
            
        - 1: only contig(C); 

            + Select this option if you only have contigs and one reference that is complete. PhaME requires users to input a complete genome as reference. This option will only report contigs in the final alignments and the tree.

        - 2: only reads(R);

            + Select this option if you only want to analyze reads file. Remember similar to option 1, a reference must be given.

        - 3: combination F+C;

            + Select this option if you have full/complete genome and contigs.

        - 4: combination F+R;

            + Select thsi option if you have full/complete genomes and only reads.

        - 5: combination C+R; 

            + Select this option if you have full/complete genomes and reads. However, it still requires at least one Full/Complete genome.

        - 6:combination F+C+R;

            + Select this option if you have full/complete, contigs, and read datasets in your analysis.

        - 7:realignment 

            + Select this option if you want to realign using a subset of genomes that have already been aligned using one of the option above. It requires editing the `working_list.txt` file.

10. *reads*
    This option is dependent on option chosen in `data`. If the analysis contains only single reads, enter 1, if paired reads enter 2, and if both are present enter 3.

11. *tree*
    The option to generate tree. If 0 is entered, no tree is generted. If 1 is entered, only FastTree is used. If 2 is entered, only RAxML is used. If 3 is entered, both FastTree and RAxMl are used to make trees.

12. *bootstrap*
  - The option is valid if 2 or 3 is selected in `tree` option. It will calculate bootstrap trees using RAxML.

13. *N*
  - Specify the number of bootstrap trees to generate if its turned ON in `bootstrap` option.

14. *PosSelect*
    The option to turn ON and select type of molecular evolution analysis to be done. Enter 0 to turn OFF molecular evolutionary analysis, 1 to use PAML to do molecular evolutionary analysis, 2 to use HyPhy, and 3 to use both of them. Turning this option ON will significantly slow the runtime. If this option is turned ON, you must provide the gff file for the corresponding reference genome.

15. *code*
    This specifies the pre-calculated parameters during genome alignments.Option 0 which is specific for bacteria uses, `Bacteria` aligns using default option with `maxmatch` for nucmer. And, option 1 which is for`Virus` sets option for nucmer alignment with `maxmatch` turned ON and `-b 200 -c 65 -d 0.12 -g 90 -l 20`.

16. *clean*
    Turning this option ON (1) will remove intermediate files.

17. *threads*
    Specify the number of threads to run analysis ON.

18. *cutoff*
    This options lets user control the genomes to include based on how much of their region was included in the alignemnt against the reference genome. Linear alignment (LA) coverage against reference - ignores SNPs from organism that have lower cutoff.


What output files are produced?
===============================
  
Summary files ( all files are found under folder `workdir/results`)

    - SNP alignment files

        + all positions (including monomorphic sites)

            * `project`_all_alignment.fna
        - all detected SNPs
            `project`_all_snp_alignment.fna
        - SNPs in CDS (coding sequence)
            `project`_cds_snp_alignment.fna
        - intergenic SNPs
            `project`_int_snp_alignment.fna

    - Newick tree files
        - bootstrap mapped Maximum Likelihood trees
            - RAxML_bipartitionsBranchLabels.`project`_all_best
            - RAxML_bipartitions.`project`_all_best 
        - bootstraps
            - RAxML_bootstrap.`project`_all_b
        - best ML tree
            - RAxML_bestTree.`project`_all
        - RAxML tree using only CDS SNPs
        - FastTree using all SNPs
            - `project`_all.fasttree
    - FastTree using SNPs from coding sequence
        `project`_cds.fasttree
    
    - Other files:
        - coordinates of gaps throughout the overall alignment
            `project`_gaps.txt
        - the size of gaps between `reference` and other genomes.
            `project`_all_gaps.txt
        - A tab delimited summary file containing information on the core genome size, total SNPs, etc.
            `project`_summaryStatistics.txt 
                - Most rows are genome name (first column), attribute name (second column), and corresponding value (third column)
                  - `Total_length` for genome size (total base pair) of the corresponding genome (first column)
                - `Gap_legnth` for total gaps in the corresponding genome (first column)
                - One row labeled `REPEAT` (first column) and `Gap_length`(second column) correspond to repeat size (third column) of reference genome.
                - `Reference used` shows the name of the reference genome used.
                - `Total gap length:` shows the length of total gaps in the alignment.
                - `Core genome length:` shows the length of genomes that were aligned.
                - `Total SNPs:` shows the length of SNPs.
                - `CDS SNPs:` shows the subset of SNPs from Total SNPs that fall within coding regions.
        - A pairwise list of all compared position with coordinates between references and samples
            `project`_comparison.txt
            `project`_stats.txt (also contains if SNPs are in coding or non-coding regions)
        - A matrix file that lists the number of SNPs present between genomes
            - all core regions
                `project`_snp_coreMatrix.txt
            - CDS only
                `project`_snp_CDSmatrix.txt
            - intergenic only
                  `project`_snp_intergenicMatrix.txt
    - Log file
        `project`.log
    
    - Error file
         `project`.error 

Directory structures	

    - `working directory`/files
         references (concatenated chromosomes)
    - `working directory`/results
         All output files
    - `working directory`/results/snps
        SNP coordinate files generated from NUCmer and bowtie
            - `g1_g2.snps`: contains pairwise snps between `g1` and `g2`. For example:
    
    .. code-block:: console

        [P1] [SUB]   [SUB]   [P2]    [BUFF]  [DIST]  [FRM]   [TAGS]
        127     T       C    127        22      127     1       1   KJ660347_1_18959    ZEBOV_2002_Ilembe_1_18958
        149     T       C    149        6       149     1       1   KJ660347_1_18959    ZEBOV_2002_Ilembe_1_18958
        155     C       A    155        6       155     1       1   KJ660347_1_18959    ZEBOV_2002_Ilembe_1_18958


  - `working directory`/results/gaps
      - Gap coordinate files generated from NUCmer and bowtie
  - `working directory`/results/stats
      - Intermediate stat files generated when parsing NUCmer and Bowtie results
        - `g1_g2.coords` is a table file that contains regions of genome `g1` and `g2` that were aligned.
        - `g_repeat_coords.txt` is a table that contains region within genome `g` that were detected as similar.
        - `g_repeat_stats.txt` contains genome size, repeat segment, and repeat length of genome `g`. For example:
        
        ::

            ZEBOV_2007_4Luebo size: 18958
            Repeats segment #:  0
            Repeats total length:   0 (0.00%)

        - `repeat_stats.txt` summary of all `g_repeat_stats.txt`.

  - `working directory`/results/temp
      - Temporary files generated
  - `working directory`/results/PSgenes
      - All gene fasta files  that contain at least 1 SNP, along with their amino acid sequences and codon alignments
  - `working directory`/results/paml
      - PAML results
  - `working directory`/results/hyphy
      - HyPhy results
  - `working directory`/results/`*_ambiguousSNPpositions.txt`
      - Positions in genomes represented as raw reads where there are ambiguous SNPs.
