import sys
from CountHexanuc import CuentaPalindromos as CountPal

try:
    GenomePath = sys.argv[1]
    CountPal(str(GenomePath))
except FileNotFoundError:
    print ("There directory path \"{}\"is incorrect".format(GenomePath))