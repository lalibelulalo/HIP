def CuentaPalindromos(GENOMES_PATH):
    import os
    import glob
    import time
    strings = time.strftime("%Y,%m,%d,%H,%M,%S")
    t = strings.split(',')
    numbers = [ int(x) for x in t ]  
    current_date = "-".join([str(numbers[0]),str(numbers[1]),str(numbers[2])])
    current_time = "".join ([str(numbers[3]),'hrs',str(numbers[4]),'mins'])
    date = "_".join ([current_date,current_time])
    output_file = "".join (['Markov_count_',str(date),'.txt'])
    output_file = str(output_file)
    print ("The output file is: {}".format(output_file))    
    
    # get the start time
    st = time.time()
    InitDir = os.getcwd()
    FilesDir = "/".join ([InitDir,GENOMES_PATH])
    output = open (output_file, 'w')
    output.write("spp,spp,palindrome,obs,markov0,markov1,markov2,markov3,genomesize,A,T,C,G,N\n")
    header = ['>']
    
    path = GENOMES_PATH
    genomes = [x for x in os.listdir(path) if x.endswith(".fasta") or x.endswith(".fna")]
    
    Nucs = ["N","A","G","C","T"]
    # CARGAMOS PALINDROMOS
    pals = [Nucs[i]+Nucs[j]+Nucs[k]+Nucs[-k]+Nucs[-j]+Nucs[-i] for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs))) for k in range(1,(len(Nucs)))]
    contador_genomes = 1
    
    try:
        for GenomeFile in genomes:
            st2 = time.time()
            print ("Archivo {} de {}: {}".format(contador_genomes,len(genomes),GenomeFile))
            GenomeFileName = GenomeFile
            GenomeFile = "".join ([GENOMES_PATH,GenomeFile])
            No_header = "".join ([GenomeFile, ".No_Header"])
            One_line = "".join ([No_header, ".One_line"])

            with open(GenomeFile) as genome_Header, open(No_header, 'w') as genome_NoHeader:
                for line in genome_Header:
                    if not any(sign in line for sign in header):
                        genome_NoHeader.write(line)
            genome_NoHeader.close()

            with open(No_header, 'r') as f, open(One_line,'w') as genome_OneLine:
                lines = f.readlines()
                mystr = ''.join([(line.strip()).lower() for line in lines])
                genome_OneLine.write(mystr)
            genome_OneLine.close()

            genome_file = open(One_line)
            genome = genome_file.read()

            # BUSCAMOS NUCLEOTIDOS EN TODO EL GENOMA Y OBTENEMOS EL TAMAÑO DEL GENOMA
            nuc_genome, nuc_list = [genome[ii:ii+1] for ii, ch in enumerate(genome)], ('a','t','c','g','n')
            NUCLEOTIDES, genome_length = {nuc:nuc_genome.count(nuc) for nuc in nuc_list}, len(nuc_genome)

            # k-meros del genoma
            dinuc_genome, trinuc_genome, tetranuc_genome, hexanuc_genome = [genome[ii:ii+2] for ii, ch in enumerate(genome)], [genome[ii:ii+3] for ii, ch in enumerate(genome)], [genome[ii:ii+4] for ii, ch in enumerate(genome)], [genome[ii:ii+6] for ii, ch in enumerate(genome)]

            #----------------------------------------------------------------
            # HACEMOS EL CONTEO PARA CADA PALINDROMO
            for pal in pals:
                pal = pal.lower()
                #print ("Palindromo: {}".format(pal))

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


                ## MODELO DE ORDEN 0:
                # Producto de frecuencias de nucleotidos del palindromo en el genoma (DENOMINADOR)
                frec_nuc = 1
                for ii, ch in enumerate(pal):
                    nucleotide = pal[ii:ii+1]
                    if not (len(nucleotide) == 1):
                        break
                    frec_nuc = NUCLEOTIDES[nucleotide]*(1/genome_length)*frec_nuc
                    # DIVISION ENTRE FRACCIONES (MODELO DE MARKOV ORDEN 1)
                markov0 = (frec_nuc)*(genome_length - 5)
                ###print ("Expected M0: {} - O/E: {}".format(markov0,(OCTANUCLEOTIDES[pal]/markov0)))

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
                markov1 = (numerador_frec_dinuc/denominador_frec_nuc)*(genome_length - 5)
                ###print ("Expected M1: {} - O/E: {}".format(markov1,(OCTANUCLEOTIDES[pal]/markov1)))

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
                markov2 = (numerador_frec_trinuc/denominador_frec_dinuc)*(genome_length - 5)
                ###print ("Expected M2: {} - O/E: {}".format(markov2,(OCTANUCLEOTIDES[pal]/markov2)))

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
                markov3 = (numerador_frec_tetranuc/denominador_frec_trinuc)*(genome_length - 5)
                ###print ("Expected M3: {} - O/E: {}".format(markov3,(OCTANUCLEOTIDES[pal]/markov3)))
                #print ("------------------------------------")     
                output.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(GenomeFileName,GenomeFileName,pal,HEXANUCLEOTIDES[pal],markov0,markov1,markov2,markov3,genome_length,NUCLEOTIDES['a'],NUCLEOTIDES['t'],NUCLEOTIDES['c'],NUCLEOTIDES['g'],NUCLEOTIDES['n']))
            # get the end time
            et1 = time.time()
            # get the execution time
            elapsed_time1 = et1 - st2
            print('Execution time: {} secs.'.format(elapsed_time1))
            print ("------------------------------------")
            contador_genomes = contador_genomes + 1
    except FileNotFoundError:
        print ("There palindrome file \"{}\" not exist.")
    output.close()
    
    for NoHeaderPath in glob.iglob(os.path.join(FilesDir, '*.No_Header')):
        os.remove(NoHeaderPath)
    One_line.close()
    for OneLine in glob.iglob(os.path.join(FilesDir, '*.One_line')):
        os.remove(OneLine)
    # get the end time
    et2 = time.time()
    # get the execution time
    elapsed_time2 = et2 - st
    print('Code execution time : {} mins.'.format(elapsed_time2/60))
    print ("***")
