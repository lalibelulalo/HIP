Introduction
#############

¿Qué es PhaME?
==============

PhaME o Phylogenetic and Molecular Evolution (PhaME) es una herramienta de análisis que realiza análisis filogenéticos y de evolución molecular.

Dada una referencia, PhaME extrae SNP de genomas completos, borradores de genomas y/o lecturas, utiliza alineación de secuencias múltiples de SNP para construir un árbol filogenético y proporciona análisis evolutivos (genes bajo selección positiva) utilizando CDS SNP.


PhaME: Bajo la capa
======================

La entrada de PhaME consiste en un conjunto de genomas en formatos fasta y/o fastq, y los archivos de anotación correspondientes en formato gff3 si los análisis posteriores incluirán regiones de codificación. A continuación se proporciona una explicación detallada paso a paso de cómo PhaME analiza los genomas.

Selección del genoma de referencia:
-----------------------------
Dado que PhaME es una herramienta basada en el genoma de referencia en la que todos los genomas y metagenomas de entrada se alinean con una referencia, el primer paso del análisis de PhaME es seleccionar un genoma de referencia. Dado un conjunto de genomas (en una carpeta enumerada bajo el parámetro "refdir" del archivo de control), el genoma de referencia se puede seleccionar mediante una de tres opciones: opción 1: se selecciona un genoma aleatorio del conjunto de genomas proporcionado; opción 2: se selecciona un genoma específico del conjunto a través de la entrada del usuario; opción 3: la distancia MinHash se calcula entre todos los genomas proporcionados (genomas completos, borradores de genomas y lecturas sin procesar) para determinar qué genoma de referencia tiene la distancia promedio más corta a todos los demás genomas. Las distancias de MinHash se calculan utilizando su implementación en BBMap [Bushnell]_.

Auto-numerización para eliminar repeticiones de genomas de referencia:
---------------------------------------------------------------
La parte de alineación del genoma de PhaME se basa en la herramienta nucmer2 para alineaciones de genomas en formato FASTA. Cada genoma incluido se alinea primero consigo mismo usando nucmer, lo que se denomina auto-numerización, y luego las regiones alineadas llamadas repeticiones se eliminan de los genomas para los análisis posteriores. El siguiente comando nucmer se usa para el paso de auto-numerización:

::

    nucmer --maxmatch --nosimplify --prefix=seq_seq ref_genomeA.fasta ref_genomeA.fasta 

::

La opción --maxmatch, que informa todas las coincidencias, se utiliza para garantizar que se informen todas las alineaciones posibles para eliminar al máximo las repeticiones.

Alineamientos de Genoma
--------------------------------
Todos los genomas que están en `refdir` se alinean primero con el genoma de referencia (ver sección 1) al que se le han eliminado las repeticiones (sección 2). Del mismo modo, los genomas incompletos o contigs, los que se enumeran en `workdir` con la extensión `.contig`, también se alinean contra el genoma de referencia usando `nucmer`. Para alinear genomas en formato fasta entre sí, se usa el siguiente comando, igual que el paso anterior para la alineación de nucmer:

::

    nucmer --maxmatch refgenome.fasta genome.fasta

::

Todas las demás opciones en las alineaciones numéricas se mantienen por defecto, algunas de las más importantes se enumeran a continuación:

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

Además, `nucmer` solo alinea `"A"` `"T"` `"C"` `"G"`, todos los demás caracteres se ignoran. Por lo tanto, si hay "N" en los genomas provistos, estas posiciones no se incluyen en la alineación.

*Nota*: Si un análisis requiere ejecutar varias iteraciones de PhaME en un mismo conjunto de datos o un subconjunto de conjuntos de datos, no es necesario realizar la alineación una y otra vez. PhaME proporciona una opción en la que puede mantener todos los posibles alineamientos por pares de genomas de `refdir` para futuros análisis. Todos los pasos mencionados en esta sección son iguales, excepto que la alineación de todos contra todos se realiza en comparación con una sola referencia. Al hacer la alineación de todos contra todos, también se puede probar el efecto en sus análisis con diferentes genomas de referencia.

Mapeo de lecturas sin procesar al genoma de referencia
-------------------------------------------
Actualmente, PhaME solo procesa lecturas cortas y sin procesar de Illumina. Si se incluyen en los análisis lecturas sin procesar, de extremo único o emparejado, se asignan al genoma de referencia utilizando bowtie2 o BWA en función de la entrada de los usuarios. Para la asignación de lecturas al genoma de referencia, se utilizan los siguientes comandos:

Primero, construye una base de datos a partir del genoma de referencia.
::

    bowtie2-build refgenome refgenome

::
o, si se eligió BWA como el alineador preferido:

::

    bwa index refgenome

::

Luego, las lecturas sin procesar se asignan al genoma de referencia utilizando uno de los siguientes comandos: Para bowtie2 y lecturas emparejadas:
::

    bowtie2 -a -x $refgenome -1 read1 -2 read2 -S paired.sam`;

::
La opción `-a` informa de todas las alineaciones posibles.

Para bowtie2 y lecturas de un solo extremo:

::

    bowtie2 -a -x $refgenome -U read -S single.sam`;

::

Para BWA y lecturas emparejadas:

::

    bwa mem refgenome read1 read2 | samtools view -ubS -| samtools sort -T tmp_folder -O BAM -o paired.bam

::

Para BWA y lecturas de un solo extremo:

::

    bwa mem refgenome read |samtools view -ubS - | samtools sort -T tmp_folder -O BAM -o single.bam

::


Filtrado de alineaciones del genoma
------------------------------
La alineación del genoma producida usando 'nucmer' se filtra usando 'delta-filter' para mantener solo las alineaciones 1 a 1, lo que permite reorganizaciones. Este paso de filtrado se produce para todas las alineaciones `nucmer`.

::

    delta-filter -1 genome.delta > genome.snpfilter

::


Llamando SNP's desde alineaciones del genoma
--------------------------------------
Las alineaciones 'nucmer' por pares se analizan para producir una tabla SNP usando 'show-snps'.

::

    show-snps -CT genome.snpfilter > genome.snps

::

Aquí, las opciones C y T especifican no informar los SNP de alineaciones ambiguas e informar la salida en un archivo delimitado por tabulaciones, respectivamente.

Informes de alineaciones numéricas
------------------------------

Cada alineación se analiza más para producir un archivo delimitado por tabuladores que tiene información sobre las regiones y el %ID de sus alineaciones.
::

    show-coords -clTr genome.snpfilter > genome.coords

::

El indicador de parámetro -clTr implica que se informarán diferentes encabezados en el informe.

::

-c          Include percent coverage information in the output
-l          Include the sequence length information in the output
-r          Sort output lines by reference IDs and coordinates
-T          Switch output to tab-delimited format

::

Llamando SNP's desde el mapeo de lectura
---------------------------------
`bcftools mpilup` se utiliza para llamar a los SNP a partir de los resultados del mapeo de lectura (archivo bam) de todos los genomas representados por lecturas sin procesar. La profundidad máxima se establece en 1000000 para las llamadas de SNP e indel y las brechas mínimas para llamar a un indel se establecen en 3. El archivo vcf de salida se usa luego para llamar a los SNP usando `bcftools call` donde la ploidía se especifica como `1` si es un genoma haploide o bacteriano, de lo contrario se llama usando el parámetro predeterminado. Además, según el parámetro especificado por el usuario en el archivo de control, los SNP se filtran aún más según el porcentaje de SNP. Estos son los fragmentos de comando que se ejecutan como parte de esto. Todos ellos dan como resultado un archivo vcf.

::

    bcftools mpileup -d 1000000 -L 1000000 -m 3 -Ov -f $refgenome $bam_output | bcftools call --ploidy 1 -cO b > $bcf_output;
    bcftools view -v snps,indels,mnps,ref,bnd,other -Ov $bcf_output | vcfutils.pl varFilter -a$min_alt_bases -d$min_depth -D$max_depth > $vcf_output`;
    bcftools filter -i '(DP4[0]+DP4[1])==0 || (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) > $snp_filter' $vcf_output > $vcf_filtered`

::


Cálculo de las alineaciones del genoma central
----------------------------------
Como primer paso para calcular el núcleo del genoma, se comprueba la cobertura lineal de todas las alineaciones con respecto a la referencia para asegurar la proporción del genoma de referencia que se usó en la alineación. Si es inferior al umbral de corte (predeterminado = 0,6) establecido en el archivo de control, ese genoma se eliminará de los análisis posteriores. Luego, el resto de las alineaciones por pares que están en formato vcf o nucmer se cotejan para calcular un genoma central. Solo se conservan las posiciones de alineación que se conservan al 100%, todas las demás posiciones se eliminan de la alineación del genoma central final. PhaME produce múltiples archivos de alineación correspondientes al genoma central, como el que tiene solo los sitios variantes (`_all_snp_alignment.fna`), tiene sitios variantes e invariantes (`all_alignment.fna`) y los que tienen SNP solo de la codificación. región (`_cds_snp_alignment.fna`). La alineación SNP de la región de codificación requiere un archivo de anotación con formato GFF.


Reconstruyendo la filogenia del genoma central
-------------------------------------
PhaME proporciona varias herramientas (RAxML [Stamatakis 2014]_, FastTree [Price 2010]_ e IQ-Tree [Nguyen 2015]_) para reconstruir la filogenia a partir de alineaciones del genoma central que tienen sitios invariantes. Si se elige la opción RAxML o FastTree, los usuarios no pueden modificar los modelos ya que están preseleccionados. Los árboles RAxML se reconstruyen usando modelos GTRGAMMAI que "GTR + Optimización de tasas de sustitución + Modelo GAMMA de heterogeneidad de tasas (se estimará el parámetro alfa)" con `I` pero con estimación para sitios invariables. FastTree solo usa el modelo GTR. IQ-TREE se ejecuta usando la opción `-m TEST` que busca el mejor modelo que se ajuste a los datos antes de reconstruir la filogenia. RAxML es la única opción disponible actualmente que también puede calcular los arranques.

Selección de genes para análisis evolutivos moleculares
-------------------------------------------------------
Para realizar análisis de selección utilizando PAML o HyPhy, se requieren alineaciones de codones de genes. Según la posición de los SNP en el genoma de referencia, si un SNP se encuentra dentro de una región de codificación y si esa región de codificación no tiene una brecha, se extraen de la alineación del genoma central. Las secuencias de nucleótidos de los genes se traducen a secuencias de proteínas, se alinean usando el programa mafft 8 y luego se traducen de nuevo a nucleótidos usando el código Perl pal2nal.pl de http://www.bork.embl.de/pal2nal/.

Análisis Evolutivo Molecular
------------------------------------

El conjunto de alineaciones de genes se utiliza para análisis evolutivos moleculares mediante PAML [Yang 2007]_ o HyPhy. Ambos paquetes pueden probar la presencia de sitios y linajes seleccionados positivamente al permitir que la relación dN/dS (ω) varíe entre sitios y linajes. El modelo de prueba adaptativa REL del sitio de sucursal para la diversificación episódica (aBSREL) en el paquete HyPhy se utiliza para detectar instancias de diversificación episódica y selección positiva. Si se selecciona PAML, se implementan los modelos anidados M1a-M2a y M7-M8. En este último caso, la prueba de razón de verosimilitud entre los modelos nulo (M1a y M8) y el modelo alternativo (M2a y M7) con un límite de significación del 5 % proporciona información sobre cómo evolucionan los genes. Los resultados de cada gen luego se resumen en una tabla que contiene información sobre si el gen está evolucionando bajo selección positiva, neutral o negativa, junto con los valores p. HyPhy se ejecuta con un modelo que busca específicamente signos de selección positiva en determinados conjuntos de genes. El análisis produce una lista de archivos JSON correspondientes a cada gen que se puede cargar en vision.hyphy.org/absrel para su posterior análisis. Optamos por proporcionar PAML como una opción; sin embargo, recomendamos usar HyPhy para proyectos grandes debido a su velocidad y resultados concisos.


References
--------------
.. [Yang 2007] Yang Z: PAML 4: phylogenetic analysis by maximum likelihood. Mol Biol Evol 2007, 24:1586-1591.
.. [Pond 2005] Pond SL, Frost SD, Muse SV: HyPhy: hypothesis testing using phylogenies. Bioinformatics 2005, 21:676-679.
.. [Kurtz 2004] Kurtz S, Phillippy A, Delcher AL, Smoot M, Shumway M, Antonescu C, Salzberg SL: Versatile and open software for comparing large genomes. Genome Biol 2004, 5:R12.
.. [Bushnell] Bushnell B: BBMap. 37.66 edition. sourceforge.net/projects/bbmap/.
.. [Stamatakis 2014] Stamatakis A: RAxML version 8: a tool for phylogenetic analysis and post- analysis of large phylogenies. Bioinformatics 2014, 30:1312-1313.
.. [Price 2010] Price MN, Dehal PS, Arkin AP: FastTree 2--approximately maximum- likelihood trees for large alignments. PLoS One 2010, 5:e9490.
.. [Nguyen 2015] Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ: IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Mol Biol Evol 2015, 32:268-274.
