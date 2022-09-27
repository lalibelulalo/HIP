def CuentaHexameros(GENOMES_PATH,MARKOV):
    from Bio import SeqIO
    import os
    import glob
    import time
    import re
    import numpy as np

    ## Obtenemos tiempo inicial
    st = time.time()
    st2 = time.time()
    ## Obtenemos la fecha para nombrar el archivo de salida
    DateTime = time.strftime("%Y,%m,%d,%H,%M,%S")
    t = DateTime.split(',')
    numbers = [ int(x) for x in t ]  
    current_date = "-".join([str(numbers[0]),str(numbers[1]),str(numbers[2])])
    current_time = "".join ([str(numbers[3]),'hrs',str(numbers[4]),'mins'])
    date = str("_".join ([current_date,current_time]))
    ## Nombre del modelo
    model = str("".join(['M',str(MARKOV)]))
    ## Nombre del archivo de salida
    output_file = str("_".join (['Markov_count',GENOMES_PATH,date,'HEXANUCS',model,'.txt']))
    
    ## Imprimimos el nombre del archivo de salida
    print ("The output file is: {}.\n".format(output_file))
    
    ## Creamos el archivo de salida
    header = str("".join(['spp\tpalindrome\tobs\tmarkov',str(MARKOV),'\tgenomesize\tA\tTh\tC\tG\tN\n']))
    output = open (output_file, 'w')
    output.write(header)
    
    ## Entramos al path proporcionado que contiene los genomas
    GenomeDir = str("".join ([GENOMES_PATH,'/']))
    ## Hacemos una lista con los genomas
    genomes = [x for x in os.listdir(GenomeDir) if x.endswith(".fasta") or x.endswith(".fna")]
    
    # CARGAMOS PALINDROMOS
    Nucs = ["N","A","G","C","T"]
    pals = [Nucs[i]+Nucs[j]+Nucs[k]+Nucs[-k]+Nucs[-j]+Nucs[-i] for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs))) for k in range(1,(len(Nucs)))]
    
    try:
        contador = 0
        for GenomeFile in genomes:
            contador += 1
            print ("Archivo {} de {}: {}".format(contador,len(genomes), GenomeFile))
            ## Obtenemos nombre del genoma en turno
            GenomeFileName = GenomeFile
            GenomeFile = "".join ([GenomeDir,GenomeFile])
            ## Extraemos las secuencias sin encabezados
            fasta_sequences = SeqIO.parse(open(GenomeFile),'fasta')
            ## Creamos un string con todo el genoma
            genome = [str(fasta.seq) for fasta in fasta_sequences]
            genome = ''.join(genome)

            # BUSCAMOS NUCLEOTIDOS EN TODO EL GENOMA Y OBTENEMOS EL TAMAÑO DEL GENOMA
            nuc_genome, nuc_list = [genome[ii:ii+1] for ii, ch in enumerate(genome)], ('A','T','C','G','N')
            NUCLEOTIDES, genome_length = {nuc:nuc_genome.count(nuc) for nuc in nuc_list}, len(nuc_genome)

            # k-meros del genoma
            dinuc_genome, trinuc_genome, tetranuc_genome, hexanuc_genome = [genome[ii:ii+2] for ii, ch in enumerate(genome)], [genome[ii:ii+3] for ii, ch in enumerate(genome)], [genome[ii:ii+4] for ii, ch in enumerate(genome)], [genome[ii:ii+6] for ii, ch in enumerate(genome)]

            #----------------------------------------------------------------
            # HACEMOS EL CONTEO PARA CADA PALINDROMO
            for pal in pals:
                pal=pal.upper()

                ## Dinucleotidos del palindromo en cuestion
                dinuc_pal = [pal[ii:ii+2] for ii, ch in enumerate(pal) if len(pal[ii:ii+2])==2 ]
                DINUCLEOTIDES = {dinuc:dinuc_genome.count(dinuc) for dinuc in dinuc_pal}

                ## Trinucleotidos del palindromo en cuestion
                trinuc_pal = [pal[ii:ii+3] for ii, ch in enumerate(pal) if len(pal[ii:ii+3])==3 ]
                TRINUCLEOTIDES = {trinuc:trinuc_genome.count(trinuc) for trinuc in trinuc_pal}

                ## Tetranucleotidos del palindromo en cuestion
                tetranuc_pal = [pal[ii:ii+4] for ii, ch in enumerate(pal) if len(pal[ii:ii+4])==4 ]
                TETRANUCLEOTIDES = {tetranuc:tetranuc_genome.count(tetranuc) for tetranuc in tetranuc_pal}

                ## Hexanucleotidos del palindromo en cuestion
                hexanuc_pal = [pal[ii:ii+6] for ii, ch in enumerate(pal) if len(pal[ii:ii+6])==6 ]
                HEXANUCLEOTIDES = {hexanuc:hexanuc_genome.count(hexanuc) for hexanuc in hexanuc_pal}
                ###print ("Observed: {}".format(OCTANUCLEOTIDES[pal]))

                if (MARKOV==0):
                    ## MODELO DE ORDEN 0:
                    # Producto de frecuencias de nucleotidos del palindromo en el genoma (DENOMINADOR)
                    frec_nuc = 1
                    for ii, ch in enumerate(pal):
                        nucleotide = pal[ii:ii+1]
                        if not (len(nucleotide) == 1):
                            break
                        frec_nuc = NUCLEOTIDES[nucleotide]*(1/genome_length)*frec_nuc
                        # DIVISION ENTRE FRACCIONES (MODELO DE MARKOV ORDEN 1)
                    markov = (frec_nuc)*(genome_length - 5)

                elif (MARKOV==1):
                    ## MODELO DE ORDEN 1:
                    # Producto de frecuencias de Dinucleotidos del palindromo en el genoma (NUMERADOR)
                    numerador_frec_dinuc = 1

                    for ii, ch in enumerate(pal):
                        dinuc = pal[ii:ii+2]
                        if not (len(dinuc) == 2):
                            break
                        numerador_frec_dinuc = DINUCLEOTIDES[dinuc]*(1/(genome_length - 1))*numerador_frec_dinuc

                    # Producto de frecuencias de nucleotidos del palindromo en el genoma (DENOMINADOR)
                    contador_extremos_nuc, denominador_frec_nuc = 1, 1

                    for ii, ch in enumerate(pal):
                        nucleotide = pal[ii+1:ii+2] #pal[ii+1:ii+1+1]
                        if (contador_extremos_nuc <= 4):
                            contador_extremos_nuc += 1
                            if not (len(nucleotide) == 1):
                                break
                            denominador_frec_nuc = NUCLEOTIDES[nucleotide]*(1/genome_length)*denominador_frec_nuc

                    # DIVISION ENTRE FRACCIONES (MODELO DE MARKOV ORDEN 1)
                    if not (denominador_frec_nuc != 0):
                        markov = 0
                    else:
                        markov = (numerador_frec_dinuc/denominador_frec_nuc)*(genome_length - 5)

                    
                elif (MARKOV==2):
                    ## MODELO DE ORDEN 2:
                    # Producto de frecuencias de trinucleotidos del palindromo en el genoma (NUMERADOR)
                    numerador_frec_trinuc = 1

                    for ii, ch in enumerate(pal):
                        trinuc = pal[ii:ii+3]
                        if not (len(trinuc) == 3):
                            break
                        numerador_frec_trinuc = TRINUCLEOTIDES[trinuc]*(1/(genome_length - 2))*numerador_frec_trinuc

                    # Producto de frecuencias de dinucleotidos del palindromo en el genoma (DENOMINADOR)
                    contador_extremos_dinuc, denominador_frec_dinuc = 1, 1

                    for ii, ch in enumerate(pal):
                        dinuc = pal[ii+1:ii+3] #pal[ii+1:ii+1+2]
                        if (contador_extremos_dinuc <= 3):
                            contador_extremos_dinuc += 1
                            if not (len(dinuc) == 2):
                                break
                            denominador_frec_dinuc = DINUCLEOTIDES[dinuc]*(1/(genome_length - 1))*denominador_frec_dinuc

                    # DIVISION ENTRE FRACCIONES (MODELO DE MARKOV ORDEN 2)
                    if not (denominador_frec_dinuc != 0):
                        markov = 0
                    else:
                        markov = (numerador_frec_trinuc/denominador_frec_dinuc)*(genome_length - 5)

                    
                elif (MARKOV==3):
                    ## MODELO DE ORDEN 3:
                    # Producto de frecuencias de tetranucleotidos del palindromo en el genoma (NUMERADOR)
                    numerador_frec_tetranuc = 1

                    for ii, ch in enumerate(pal):
                        tetranuc = pal[ii:ii+4]
                        if not (len(tetranuc) == 4):
                            break
                        numerador_frec_tetranuc = TETRANUCLEOTIDES[tetranuc]*(1/(genome_length - 3))*numerador_frec_tetranuc

                    # Producto de frecuencias de trinucleotidos del palindromo en el genoma (DENOMINADOR)
                    contador_extremos_trinuc, denominador_frec_trinuc = 1, 1

                    for ii, ch in enumerate(pal):
                        trinuc = pal[ii+1:ii+4] #pal[ii+1:ii+1+3]
                        if (contador_extremos_trinuc <= 2):
                            contador_extremos_trinuc += 1
                            if not (len(trinuc) == 3):
                                break
                            denominador_frec_trinuc = TRINUCLEOTIDES[trinuc]*(1/(genome_length - 2))*denominador_frec_trinuc

                    # DIVISION ENTRE FRACCIONES (MODELO DE MARKOV ORDEN 3)
                    if not (denominador_frec_trinuc != 0):
                        markov = 0
                    else:
                        markov = (numerador_frec_tetranuc/denominador_frec_trinuc)*(genome_length - 5)
                
                else:
                    print("EL GRADO DEL MODELO ES INCORRECTO\n")
                    
                output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(GenomeFileName,pal,HEXANUCLEOTIDES[pal],markov,genome_length,NUCLEOTIDES['A'],NUCLEOTIDES['T'],NUCLEOTIDES['C'],NUCLEOTIDES['G'],NUCLEOTIDES['N']))
            
            # get the end time
            et = time.time()
            # get the execution time
            elapsed_time2 = et - st2
            st2 = time.time()
            print('Execution time: {} seconds.'.format(elapsed_time2))
            print ("------------------------------------")

    except FileNotFoundError:
        print ("POSIBLE ERROR EN PALINDROMOS.")  
    output.close()
    
    elapsed_time = et - st
    print('*** Code execution time : {} mins. ***'.format(elapsed_time/60))
