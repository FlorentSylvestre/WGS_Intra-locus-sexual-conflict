library(data.table)
###Usage:
#Rscript script infile.012 popmap target size

args = commandArgs(trailingOnly = T)
infile = args[1]
popmap = args[2]
target = args[3]
size = as.numeric(args[4])
output = "98_metrics/intersex_snp_count_differences.txt"


dat <- data.frame(fread(infile, h = F)[,-1])
pos <- fread(paste0(infile,".pos"),h = F)
colnames(dat) <- paste(pos$V1,pos$V2, sep ="_")
indv <- fread(paste0(infile,".indv"),h = F)
popmap <- fread(popmap,header = F)
popmap <- merge(indv,popmap,"V1", sort = F)

min_chr <- tapply(pos$V1, pos$V1, function(x){min(which(pos$V1 == unique(x)))})
max_chr <- tapply(pos$V1, pos$V1, function(x){max(which(pos$V1 == unique(x)))})

dat_per_sex <- split(dat,f = popmap$V2)
var_count_sex <- lapply(dat_per_sex,
                        function(x){
                          apply(x,
                                2,
                                function(y){genL = unique(y[y!= (-1)])
                                  if(length(genL) ==1 & (genL[1] %in% c(0,2))){return(0)} else{return(1)}

                                  })
                          })
rm(dat,dat_per_sex)

global = c(sum(var_count_sex[["F"]]),sum(var_count_sex[["M"]]))

list_target <- fread(target)
win = list()
pos <- paste(pos$V1,pos$V2,sep="_")
targ <- paste(list_target$V1,list_target$V2,sep ="_")

win= list()
for(i in 1:nrow(list_target)){
  target <- which(pos == targ[i])

  inf <- max(min_chr[names(min_chr) == list_target$V1[i]], target - size)
  sup <- min(max_chr[names(max_chr) == list_target$V1[i]], target + size)

  FC = var_count_sex[["F"]][inf:sup]
  MC = var_count_sex[["M"]][inf:sup]
  Wrange = range(as.numeric(unlist(strsplit(names(FC),"_"))[c(F,T)]))

  Wmidpos = sum(Wrange)/2
  Wsize = Wrange[2]-Wrange[1]
  Nf = sum(FC)
  Nm = sum(MC)
  p_fish = fisher.test(cbind(global,c(Nf,Nm)))$p.value
  win[[i]] = c(Wmidpos,Wsize,Nf,Nm,p_fish)
 if(i%%1000 == 0){print(i)}
}

res <- cbind(list_target,do.call(rbind,win))
warnings()
colnames(res) <- c("CHROM","POS", "midpos","Wsize","Nf","Nm","P_fish")


fwrite(res,output, row.names = FALSE,sep = "\t")

