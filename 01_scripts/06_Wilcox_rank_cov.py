#!/usr/bin/env python3
"""For each SNPs, report wilcoxon rank sum test p-value comparing male and female coverage
Usage:
        <program> input_vcf popmap output

"""


# Modules
from collections import defaultdict
import sys
import gzip
from scipy.stats import mannwhitneyu  as wilkox

# Parsing user input
try:
        input_vcf = sys.argv[1]
        input_map = sys.argv[2]
        output = sys.argv[3]
except:
        print(__doc__)
        sys.exit(1)

#Functions
def myopen(_file, mode = "rt"):
        if _file.endswith(".gz"):
                return gzip.open(_file, mode)
        else:
                return open(_file, mode)

def parsing_map(map_path, ind_list):

        mapfile = open(map_path,"r")
        map = defaultdict(list)

        for line in mapfile:

                infos = line.strip().split("\t")
                map[infos[1]].append(ind_list.index(infos[0]))

        return map

#Main script
with myopen(input_vcf) as inputfile:
        with myopen(output,"wt") as outputfile:
                outputfile.write("CHROM\tPOS\tstat\tpval\n")
                for line in inputfile:
                        l = line.strip().split("\t")

                        if line.startswith("#"):
                                if line.startswith("#CHROM"):
                                        popmap = parsing_map(input_map, l[9:])


                        else:

                                infos_cov = l[8].split(":").index("DP")

                                genotype = []
                                for sex in ["M","F"]:
                                        genoS = [l[9:][x] for x in popmap[sex] if l[9:][x].split(":")[0] != "./."]
                                        genotype.append(genoS)

                                stat,pval = wilkox([int(x.split(":")[infos_cov]) for x in genotype[0]], [int(x.split(":")[infos_cov]) for x in genotype[1]], alternative = "two-sided")
                                new_line = "\t".join(l[0:2])
                                new_line += "\t" + str(stat) + "\t" + str(pval) + "\n"
                                outputfile.write(new_line)


