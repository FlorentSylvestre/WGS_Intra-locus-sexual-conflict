#!bin/python3

"""
Script that perform intersex SNP-by-SNP fst calculation in sliding windows of size Size and overlap Overlap
It repports the ratio of the number of significant SNPs in real dataset over this number in permuted dataset
Averaged over N permutation, as well as the associated p-value based on permutation test P(Nperm >= Ntrue)
Its based on a vcf of which the first M sampels represent males, and other females.
!!!!!!!!!!!!Works with numpy v 1.23.4, not sure about others!!!!!!!!!!!!
Usage:
python3 permut_fst_windows.py <Vcf> <Size> <Overlap> <N> <M>

"""


from collections import defaultdict
import numpy as np
import sys
import gzip
from scipy import stats
import numpy.ma as ma

#Parsing user input:

args = sys.argv

Vcf = args[1]
Size = int(args[2])
Overlap = int(args[3])
N = int(args[4])
M = int(args[5])

if Size/Overlap != Size // Overlap:
    print("Size must be a multiple of Overlap")
    exit(0)

##Fixed varaibles
Output = "08_Fst_perm/permuted_Fst_win"


#Functions:
#Global interaction with files:
def myopen(_file, mode = "rt"):
    if _file.endswith(".gz"):
        print(mode)
        return gzip.open(_file, mode)

    else:
        return open(_file, mode)

def Permuts(N, Total):
    rng = np.random.default_rng()

    return rng.permuted(np.tile(range(Total),(N,1)), axis = 1)
  
  def Conv_geno(genos):
    genos[genos == "./."] = -9
    genos[genos == "0/0"] = 1
    genos[genos == "0/1"] = 0.5
    genos[genos == "1/1"] = 0
    return genos.astype("float")


#Stats calculations:

def Fst(M,F, Nm, Nf):
    #pm = ((M == "0/0").sum(axis = 1) + 0.5 * (M == "0/1").sum(axis = 1))/Nm
    #pf = ((F == "0/0").sum(axis = 1) + 0.5 * (F == "0/1").sum(axis = 1))/Nf
    pm = M.sum(axis = 1)/Nm
    pf = F.sum(axis = 1)/Nf
    return ((pm-pf)**2)/(4*((pm+pf)/2)*(1-((pm+pf)/2)))


def Scaling_Fst_Pvalue_factor(Nm, Nf):
    return 1 / (2 * (2 / ((1 / Nm) + (1 / Nf))))


def P_value_Fst(geno_M, geno_F):
    Nm = geno_M.count(axis = 1)
    Nf = geno_F.count(axis = 1)
    return stats.chi2.sf(Fst(geno_M, geno_F, Nm, Nf )/ Scaling_Fst_Pvalue_factor(2*Nm, 2*Nf), 1)


##Core script
#Skipping header and defining windows
with myopen(Vcf) as In:
    for line in In:
        if line.startswith("#CHROM"):
            Total = len(line.split("\t")[9:])
            Windows =  defaultdict(lambda: defaultdict(list))
            continue
        if line.startswith("#"):
            continue

        l = line.strip().split("\t")
        Genotypes = [x.split(":")[0] for x in l[9:]]
        chrom, pos = l[0:2]
        Windows[chrom][int(pos)//Overlap].append(Conv_geno(np.array(Genotypes)))
#Analyses:

Out = open(Output, "w")
step = Size // Overlap #Number of windows to aggregate
for chrom in Windows:
    maxwin = max(Windows[chrom].keys())
    for win in range(maxwin+1):
        print(f"running {chrom}: {win}")
        #Real dataset
        
                subdat = []
        [subdat.extend(Windows[chrom][x]) for x in range(win,win + step) if x in Windows[chrom]]
        NSNPs = len(subdat)

        subdat = ma.masked_values(subdat, -9)

        if NSNPs == 0:
            continue

        Count_signif_SNPs = P_value_Fst(subdat[:,:M], subdat[:,M:])
        Count_signif_SNPs = [x for x in Count_signif_SNPs if x <=0.01]
        Count_signif_SNPs = len(Count_signif_SNPs)

##Permutations i
       ##Permutations
        rand_index = Permuts(N,Total)
        data_permute = [subdat[:,x] for x in rand_index]
        data_permute = ma.array(data_permute).reshape((NSNPs * N, Total))

        p_values_perm = P_value_Fst(data_permute[:,:M], data_permute[:,M:])

        permend = time.time()

        Res_count_perm = [0] *N
        for x in range(0,N):
            count_perm = [x for x in p_values_perm[x*NSNPs:(x+1)*NSNPs] if x <= 0.01]
            Res_count_perm[x] = len(count_perm)

        ratios = [Count_signif_SNPs / x for x in Res_count_perm if x != 0]
        mean_Ratio = sum(ratios) / len(ratios) if len(ratios) > 0 else "NA"
        p_value = [x for x in Res_count_perm if x >= Count_signif_SNPs]
        Out.write(str(chrom)  + "\t" + "\t".join([str(win*Overlap), str(min(maxwin*Overlap, (win+ Size//Overlap)*Overlap))]) + "\t"  +str(NSNPs) + "\t" + str(Count_signif_SNPs) + "\t" + str(len(ratios)) + "\t" + str(mean_Ratio) + "\t" + str(len(p_value)/N) + "\n")

                                                                           
  
