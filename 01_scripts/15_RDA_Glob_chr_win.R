#!/usr/bin/env Rscript

#Script that perform RDA analysis using sex as only constraining variable
#Perform RDa at :
#1) global scale
#2) per chr scale
#3) per sliding windows of size N

#Usage : Rscript RDA_win.R Infile popmap windows_size 
#with:
# Infile : 012 file as outputeed by vcftools.  ex : "filtered_data.012". will read 012,012.indv and 012.pos files
# popmap : True sex affiliation for all individual. 2 tab separated column : first == sample code, second = M or F (Male or Female)
# windows_size = windows size in bp



##parsing user arguments
args=commandArgs(trailingOnly=TRUE)
Infile = args[1]
popmap = args[2]
size = as.numeric(args[3])


output = "09_RDA/OUTPUT_RDA"
step = 0.2
nthread = 10



###packages needed
library(data.table)
library(vegan)

###Functions
define_windows <- function(pos, size){
    #Define windows with 50% overlap
    midpos <- size/2 + size * (pos %/%size)
    stand_pos <- (pos - size/2)
    midpos_overlap <- size + size * (stand_pos%/%size)
    midpos_overlap[midpos_overlap == 0] <- NA
    return(cbind(midpos,midpos_overlap))}



###Parsing popmap:
print("parsing map")
maps = read.table(popmap, h = F)
colnames(maps) <- c("Indv", "Real_Sex")

##ordering it according to order in 012 file:

indv <- read.table(paste0(Infile,".indv"))
colnames(indv) <- "Indv"
maps <- base::merge(indv, maps,by = "Indv", sort = F)
print(maps)

###reading posfile and preparing windows:
print("preparing windows")
pos = fread(paste0(Infile, ".pos"))
colnames(pos) = c("CHROM", "POS")
pos$CHROM = factor(pos$CHROM, levels = unique(pos$CHROM))

###reading input file
dat <- as.data.frame(fread(Infile, header = F)[,-1])
print(dim(dat))
print(dim(pos))
colnames(dat) <- paste(pos$CHROM,pos$POS, sep = "_")
rownames(dat) <- indv$Indv

###Initializing output list
Res_score = list()
Res_Signif = list()


###Global RDA
print("Global RDA")
RDA_global = rda(dat ~ maps$Real_Sex, scale = F)
type = "Global"

Res_score[[type]] = cbind(type,
    pos$CHROM,
    pos$POS,
   "NA",
   "NA",
    scores(RDA_global, choices = c(1), display = "species"))
print(Res_score[[type]])
colnames(Res_score[[type]]) = c("Scale","chr","Pos","Start","End","Contribution")
pval_rda = anova.cca(RDA_global,permutations = 1000, parallel = nthread)$'Pr(>F)'[1]
R2_rda = RsquareAdj(RDA_global)
print("signif")
Res_Signif[[type]] = cbind(type,"NA",
    pval_rda,
    R2_rda$r.squared,
    R2_rda$adj.r.squared)
colnames(Res_Signif[[type]]) = c("Scale","Pos","P-value","R2","R2-adj")


##Per chr RDA
CHR_Obs = list()
type = "PER_CHR"

for(chr in unique(pos$CHROM)){

    CHR_Obs[[paste(chr)]] = rda(dat[,pos$CHROM == chr] ~ maps[,2])

    Res_score[[chr]] = cbind(type,
        chr,
        pos$POS[as.character(pos$CHROM) == chr],
    "NA",
    "NA",
        scores(CHR_Obs[[chr]], choices = c(1), display = "species"))
            pval_rda = anova.cca(CHR_Obs[[chr]],permutations =1000, parallel = nthread)$'Pr(>F)'[1]

    R2_rda = RsquareAdj(CHR_Obs[[chr]])

    Res_Signif[[chr]] = cbind(type,chr,
    pval_rda,
    R2_rda$'r.squared',
    R2_rda$'adj.r.squared')

}
print("CHR Done")

###Per Windows RDA
print("yeahboy")
RDA_win = list()
type = "PER_WIN"


for(chr in unique(pos$CHROM)){
    start = 0
    end = size
    while(start <= max(pos$POS[pos$CHROM == chr])){

        subdat = dat[,pos$CHROM == chr & pos$POS >= start & pos$POS <= end]
        if(ncol(subdat) !=  0){
        if(start%%size *1000 == 0){print(paste(chr,start))}
            midpos = paste(chr,(start + end)/2, sep = "_")
            RDA_win[[midpos]] = rda(subdat ~ maps[,2])

            Res_score[[midpos]] = cbind(type,
                chr,
                pos$POS[pos$CHROM == chr & pos$POS >= start & pos$POS <= end],
            start,
            end,
                scores(RDA_win[[midpos]], choices = c(1), display = "species"))

            pval_rda = anova.cca(RDA_win[[midpos]],permutations = 20000, parallel = nthread)$'Pr(>F)'[1]
            R2_rda = RsquareAdj(RDA_win[[midpos]])
            Res_Signif[[midpos]] = cbind(type,midpos,
               pval_rda,
               R2_rda$'r.squared',
               R2_rda$'adj.r.squared')


    }
        start = start + step * size
        end = end + step * size


    }
    }
print("all done")

###################HANDLE OUTPUT
fwrite(do.call(rbind,Res_score), paste(output,"contrib",sep = "."), row.names = F, sep ="\t", quote = F)
fwrite(do.call(rbind,Res_Signif), paste(output,"pval",sep = "."), row.names = F, sep ="\t", quote = F)
