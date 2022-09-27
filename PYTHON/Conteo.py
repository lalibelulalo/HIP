## use:
## python Conteo.py genomes_dir KMERS MODEL
## Where:
## KMERS: 
##	A. Octamers
##	B. HIP OCTAMERS
##	C. HEXAMERS
## MODEL: 0,1,2,3
## python Conteo.py all_pico_2022_fna B 3
import sys
from CountOctanucs import CuentaOctameros as CountOct
from CountHexanucs import CuentaHexameros as CountHex

try:
    GenomePath = sys.argv[1]
    Pals = sys.argv[2]
    Markov = sys.argv[3]
    if Pals == "A" or Pals == "B":
        CountOct(str(GenomePath),str(Pals),int(Markov))
    elif Pals == "C":
        CountHex(str(GenomePath),int(Markov))
    else:
        print("Option Error")
except FileNotFoundError:
    print ("There directory path \"{}\"is incorrect".format(GenomePath))
