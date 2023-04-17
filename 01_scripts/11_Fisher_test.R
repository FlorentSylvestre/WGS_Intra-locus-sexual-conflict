#perform per snps intersex fisher test on allele count to detect signifcantly differentiated snps between male and female
# need a genotype file format 0 1 2 (from vcftools for exemple) and a popmap
#popmap is a tab delimieted file indicating samples of first collumn and sex (M or F) on second
#Usage : Rscript Fishert_test.R genofile popmap, no header


library(data.table)

genofile <- commandArgs(trailingOnly = TRUE)[1]
popfile <- commandArgs(trailingOnly = TRUE)[2]

genotype <- data.frame(fread(genofile, stringsAsFactors = FALSE,header = FALSE))[,-1]
genopos <- data.frame(fread(paste(genofile,".pos",sep=""), stringsAsFactors = FALSE,header = FALSE))
genoind <- data.frame(fread(paste(genofile,".indv",sep=""), stringsAsFactors = FALSE,header = FALSE))
popmap <- data.frame(fread(popfile, stringsAsFactors = FALSE,header = FALSE))

fisher_function <- function(genotype, sex){
  genocount <- matrix(data = rep(0,4),ncol=2)
  colnames(genocount) <- c("M","F")
  rownames(genocount) <- c("A","R")

  for(i in 1:length(genotype)){
    if(genotype[i] == -1){
      next()
    } else if(genotype[i] == 0){
      genocount["R", sex[i] ] <-  genocount["R", sex[i] ] + 2

    }else if(genotype[i] == 2){
      genocount["A", sex[i] ] <-  genocount["A", sex[i] ] + 2

    }else if(genotype[i] == 1){
      genocount["A", sex[i] ] <-  genocount["A", sex[i] ] + 1
      genocount["R", sex[i] ] <-  genocount["R", sex[i] ] + 1

    }
  }
  if(sum(genocount[1,]) == 0) { return(NA)}else{
      return(fisher.test(genocount)$p.value)}

}

####matchimg ind position and sex from the map

genoind$Sex <- merge(genoind,popmap,sort = FALSE)[,2]

####calculing p-value for each loci
genopos$p_value <- unlist(apply(genotype, 2 , fisher_function,sex = genoind$Sex))
genopos$q_value <- p.adjust(genopos$p_value, method = "BH")

write.table(genopos, file = "07_fisher_test/fishser_test.pval",sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE)
