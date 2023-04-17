#!/usr/bin/env python3
"""Calculates per snps intersex cov and heterozogoty biais and repport alt hommozygotie
Usage:
        <program> input_vcf popmap

"""


# Modules
from collections import defaultdict
import sys
import gzip
from statistics import median

# Parsing user input
try:
        input_vcf = sys.argv[1]
        input_map = sys.argv[2]
except:
        print(__doc__)
        sys.exit(1)
        
#Fixed Path:
output = "98_metrics/Intersex_covhet"

#Functions
def myopen(_file, mode = "rt"):
        if _file.endswith(".gz"):
                return gzip.open(_file, mode)
        else:
                return open(_file, mode)


def Parse_hom(geno_list):
        geno = [x.split(":")[0] for x in geno_list]

        hom = geno.count("1/1")
        return hom


def Parse_het(geno_list):
        geno = [x.split(":")[0] for x in geno_list]
        het = geno.count("0/1")

        return het


def Parse_cov(geno_list, poscov):
        geno = [int(x.split(":")[poscov]) for x in geno_list]
        cov = median(geno)

        return cov


def parsing_map(map_path, ind_list):

        mapfile = open(map_path,"r")
        map = defaultdict(list)

        for line in mapfile:

                infos = line.strip().split("\t")
                map[infos[1]].append(ind_list.index(infos[0]))

        return map

def conv_hom(sumstats):
        Nhom = sum(sumstats["Hom"])
        Nhet = sum(sumstats["Het"])
        N = sum(sumstats["N"])

        #print((2*Nhom+ Nhet)/(2*N))
        return sumstats["Hom"] if (2*Nhom+ Nhet)/(2*N) <= 0.5 else [sumstats["N"][x] - sumstats["Hom"][x] - sumstats["Het"][x] for x in [0,1]]

#Main script



### parsing vcf
with myopen(input_vcf) as inputfile:
        with open(output,"wt") as outputfile:
                outputfile.write("CHROM\tPOS\tSEX\tN\tHet\thom\tMedian_cov\n")
                for line in inputfile:
                        l = line.strip().split("\t")

                        if line.startswith("#"):
                                if line.startswith("#CHROM"):
                                        popmap = parsing_map(input_map, l[9:])


                        else:

                                infos_cov = l[8].split(":").index("DP")
                                sumstats = defaultdict(list)

                                genotype = []
                                for sex in ["M","F"]:
                                        genoS = [l[9:][x] for x in popmap[sex] if l[9:][x].split(":")[0] != "./."]
                                        genotype.append(genoS)


                                for i in genotype:
                                        sumstats["Het"].append(Parse_het(i))
                                        sumstats["Hom"].append(Parse_hom(i))
                                        sumstats["cov"].append(Parse_cov(i, infos_cov))
                                        sumstats["N"].append(len(i))

                                hom = conv_hom(sumstats)
                                #print(sumstats)
                                #print(hom)
                                #int("A")

                                for x,y in enumerate(["M","F"]):
                                        new_line = "\t".join(l[0:2])
                                        new_line += "\t" + y + "\t" + str(sumstats["N"][x]) + "\t" + str(sumstats["Het"][x]) + "\t" + str(hom[x]) + "\t" + str(sumstats["cov"][x]) +"\n"
                                        outputfile.write(new_line)

