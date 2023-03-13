#!/usr/bin/env python3

"""input vcf using the most frequent genotype globally
Usage <script> vcf_in vcf_out
"""
import re
import sys
import gzip
from collections import defaultdict




def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def impute_simple(Lgeno):
    Main_geno= ""
    Count = 0
    for x in  set([y for y in Lgeno if y != "./."]):
        if Lgeno.count(x) > Count:
            Count = Lgeno.count(x)
            Main_geno = x

    return [Main_geno if x == "./." else x for x in Lgeno]



# Parsing user input
try:
    vcf_path = sys.argv[1]
    vcf_out  = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

### parsing first 4 column
with myopen(vcf_path) as inputfile:
    with myopen(vcf_out,"wt") as outputfile:

        for line in inputfile:
            l = line.strip().split("\t")

            if line.startswith("#"):
                    outputfile.write(line)
            else:
                start_line = "\t".join(l[:9])

                geno = [x.split(":") for x in l[9:]]
                geno_imputed = impute_simple([x[0] for x in geno])

                Ngeno = [":".join([geno_imputed[x]] + geno[x][1:]) for x in range(len(geno))]

                outputfile.write(start_line + "\t" + "\t".join(Ngeno) + "\n")
