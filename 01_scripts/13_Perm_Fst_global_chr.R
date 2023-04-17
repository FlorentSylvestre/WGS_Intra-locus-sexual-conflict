#!/usr/bin/env Rscript
#Package
library(tidyverse)
library(data.table)
library(purrr)
library(magrittr)

#Functions

nei_fst_func <- function(pm,pf){((pm-pf)^2)/(4*((pm+pf)/2)*(1-((pm+pf)/2)))}
#Parameter of the chi-square distribution given nm and nf
exp_chisq <- function(nm,nf){return(1 / ( 2 * ( 2 / ((1/nm)+(1/nf)) ) ))}
# Small function that calculate the quantile of the expected null distribution in which a SNP fall, given nm and nf and the observed allele frequencies. Usefull to compare SNP with different number of males and females
quantile_from_expdistrib<-function(pm,pf,nm,nf){
  exp_snp<-exp_chisq(nm,nf)
  quantilesnp <- pchisq(nei_fst_func(pm,pf)/exp_snp,df=1)
  return(quantilesnp)
}


output_Fst_Pval<-function(data_observed){

  #calculate FST from allele frequencies
  fst.obs<-nei_fst_func(data_observed$male_freq,data_observed$female_freq)

  #calculate in which quantile each SNP falls respective to its own theoretical distribution
  quantile_observed<-quantile_from_expdistrib(data_observed$male_freq,data_observed$female_freq,data_observed$n_males_allele_covered,data_observed$n_females_allele_covered)


  return(cbind(fst.obs,quantile_observed))
}

Count_perc_quantile<-function(quantile_observed){

  #calculate FST from allele frequencies
  #fill a vector with the number of obseerved SNPs that fall in each quantile
  quantiles.perc<-c()
  for (i in 1:100){
    quantiles.perc[i] <- length(which(quantile_observed<i/100. & quantile_observed>=(i-1)/100.))
  }

  return(quantiles.perc)
}

options(scipen=999)

args = commandArgs(trailingOnly = T)
print(args)
file_observed_data <-  args[1]

N <- as.numeric(args[3])
size = as.numeric(args[4])

#Fixed variables
suffix = "96_permutations/R_"
output = "08_Fst_perm/"


file <- fread(file_observed_data, h = T)
file$scaf <- factor(file$scaf, levels = unique(file$scaf))
file <- cbind(file,output_Fst_Pval(file))


print("Global Fst distrib for real data")
global_perc_quantile <- Count_perc_quantile(file$quantile_observed)

print("Per CHR distrib for real data")
per_chr_perc_quantile <- do.call(rbind, tapply(file$quantile_observed, file$scaf, Count_perc_quantile))



####reading permutations

res_matrix_global <- list()
res_matrix_per_chr <- list()

for(j in 1:N){
  print(paste("Perm",j,sep =" "))
  name = paste(suffix,j,sep = "")
  print(name)
  tmp_file = fread(name, h = T)
  tmp_file$scaf <- file$scaf
  tmp_file <-cbind(do.call(rbind,tapply(tmp_file$pos, tmp_file$scaf, define_windows, size = size, win = F)),tmp_file)
  tmp_file <- cbind(tmp_file,output_Fst_Pval(tmp_file))
  file <- cbind(file, tmp_file[,c("fst.obs", "quantile_observed")])

  print("Global")
  res_matrix_global[[j]] <- Count_perc_quantile(tmp_file$quantile_observed)
  print("per CHR")
  res_matrix_per_chr[[j]] <-  do.call(rbind, tapply(tmp_file$quantile_observed, tmp_file$scaf, Count_perc_quantile))
  print(res_matrix_per_chr[[j]])

  }


####Outputing files
##Global file of Fst/Pval calculation:
print("Outputing global results")
write.table(file,file = paste0(output,prefix,"_fst_res_per_snps.txt", sep =""),quote = F,sep = "\t",na = "NA", row.names = F,col.names = T)

##Global file count:
Global_res = cbind("Obs" = global_perc_quantile, do.call(cbind,res_matrix_global))

colnames(Global_res)[-1] = paste0("R",1:N)
write.table(Global_res,file = paste0(output,p refix, "_global_count_percentile.txt", sep =""),quote = F,sep = "\t",na = "NA", row.names = F,col.names = T)

##per chr count:
print("Outputing per chr result")
chr_names <- row.names(per_chr_perc_quantile)

res_per_chr <- cbind(chr_names, "Obs", per_chr_perc_quantile)
colnames(res_per_chr) <- c("CHROM", "Iter", 1:100)

for(i in 1:N){
  res_per_chr <- rbind(res_per_chr, cbind(chr_names, i,res_matrix_per_chr[[i]]))

}

write.table(res_per_chr,file = paste0(output, prefix, "_per_chr_count_percentile.txt", sep =""),quote = F,sep = "\t",na = "NA", row.names = F, col.names = T)


