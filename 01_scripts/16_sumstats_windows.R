winsize = 250000
step = 50000
thresh = 3

#function
count_Q_Fst<-function(X,Q){ apply(X,2,function(x,Q){return(sum(x>=Q,na.rm = TRUE))},Q = Q)}

chrconv <- data.frame("ARABE" = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21),
                      "ROMAIN" = c("chrI","chrII","chrIII","chrIV","chrV",
                                   "chrVI", "chrVII","chrVIII", "chrIX","chrX","chrXI","chrXII",
                                   "chrXIII","chrXIV","chrXV","chrXVI","chrXVII",
                                   "chrXVIII", "chrXX", "chrXXI"))


#Packages
library(data.table)
library(tidyr)
library(magrittr)
library(dplyr)
library(purrr)
options(scipen = 999)

#Gathering  results
RDA_contrib <- fread("09_RDA/OUTPUT_RDA.contrib")
RDA_pval <- fread("09_RDA/OUTPUT_RDA.pval")

Fst <- fread("08_Fst_perm/permuted_Fst_win", header = F)
colnames(Fst) <- c("chr","Start","End","NSNPs","Nsignif","Nratio","Mean_ratio","Pval_ratio")

theta <- fread("10_thetaNo_max_depth_250000.thetaswindow.pestPG")

Crossover <- fread("02_maps/genetic_maps.txt", skip = "set")


##parsing RDA
NSNPS <- RDA_contrib %>%
  filter(Scale == "PER_WIN") %>%
  group_by(.,across(all_of(c("chr","Start", "End")))) %>%
  summarise(nrow = length(Contribution))


Res <- data.table(NSNPS,
                  "RDA_P_val" = RDA_pval %>%
                    filter(Scale == "PER_WIN") %>% pull('P-value'),
                  "RDA_R2" = RDA_pval %>%
                    filter(Scale == "PER_WIN") %>% pull('R2'))

colnames(Res)[4] <- "NSNPS"

####Fst per win
for(chr in unique(Res$chr)){
        maxpos= max(Fst$End[Fst$chr == chr])

        Res$End[Res$chr == chr & Res$End > maxpos] <- maxpos

}
Res <- merge(Res, Fst[,-"NSNPs"])


#Theta
for(i in 1:nrow(theta)){
  theta$Chr[i]<-  chrconv$ARABE[chrconv$ROMAIN == theta$Chr[i]]
}
theta$winame <- paste(theta$Chr, theta$WinCenter, sep ="_")

Res$WinCenter <- (Res$Start + Res$End)/2
Res$winame <- paste(Res$chr, Res$WinCenter, sep ="_")
Res <- theta %>% mutate("pi" = tP/nSites) %>%
  select(nSites,Tajima,pi,winame) %>%
  inner_join(Res,.)

rm(theta)

########Recombination rates
##data
Crossover <- split(Crossover, Crossover$map)
res = list()
count = 0
for (chr in names(Crossover)) {

  if(chr != "chrUn"){maxpos = max(Fst$End[Fst$chr == chrconv$ARABE[chrconv$ROMAIN == chr]])} else{maxpos = 1000000}
  N <- max(Crossover[[chr]]$mkr)
  start = 0
  end = start + winsize

  res[[chr]] <- list()
  while(start < N){

    snps_list <- which(Crossover[[chr]]$mkr >= start  & Crossover[[chr]]$mkr <= end)
    subdat <- Crossover[[chr]][snps_list,]
    if(length(snps_list) >0){
      res[[chr]][[paste((start+end)/2)]] <- c(chr,start,min(end,maxpos), nrow(subdat),mean(subdat$loess, na.rm = T), sd(subdat$loess,na.rm = T))
    }else(count = count +1)
    start = start + step
    end = end + step

  }
  res[[chr]] <- do.call(rbind, res[[chr]])
}

res <- do.call(rbind,res) %>% as.data.frame()
res$V2 <- as.numeric(as.character(res$V2))
res$V3 <- as.numeric(as.character(res$V3))
res$V4 <- as.numeric(as.character(res$V4))
res$V5 <- as.numeric(as.character(res$V5))
res$V6 <- as.numeric(as.character(res$V6))
colnames(res) <- c("chr","Start","End","Nmark","r2","sd")
res %<>%
  filter(chr != "chrUn") %>%
  mutate("chr" = factor(res$chr[res$chr != "chrUn"], levels = chrconv$ROMAIN))

for(i in 1:nrow(res)){
  res$Chr[i]<-  chrconv$ARABE[chrconv$ROMAIN == res$chr[i]]
}

res$winame <- paste(res$Chr, (res$Start + res$End)/2, sep ="_")
res %<>%
  filter(!is.na(r2))
res$r2[res$r2 < 0] <- 0

Res <- res %>% select(winame,r2,Nmark) %>%
  inner_join(Res,.)

#Roh group
roh_group <- rep(NA, nrow(Res))
for(i in 1:nrow(Res)){
  if(Res$r2[i] <= 1) {
    roh_group[i] <- "low"
  }else if(Res$r2[i] <=5){
    roh_group[i] <- "medium"
  }else{
    roh_group[i] <- "high"
  }
}

Res$roh_group <- factor(roh_group, levels = c("low", "medium","high"))

write.table(x = Res,
            file = "11_Res_permuts/Res_permut.txt",
            row.names =  FALSE,
            sep ="\t",
            quote = FALSE)

####High Tajima regions:

Balancing_selection_windows <- Res %>%
        filter( Tajima >=quantile(Res$Tajima, 0.99)) %>%
        filter( RDA_P_val <= 0.05 | Pval_ratio <= 0.05) %>%
         select(chr,Start,End)


#merging adjascent windows
list_windows = list()
for(i in 1:nrow(Balancing_selection_windows)){
  if(i == 1){
    chr = Balancing_selection_windows$chr[i]
    start = Balancing_selection_windows$Start[i]
    end = Balancing_selection_windows$End[i]
  }else{
    if(Balancing_selection_windows$chr[i] == chr & Balancing_selection_windows$Start[i] <= end){
      end <- Balancing_selection_windows$End[i]
      if(i == nrow(Balancing_selection_windows)){list_windows[[i]] <- c(chr,start,end) }
    }else{
      list_windows[[i]] <- c(chr,start,end)
      chr = Balancing_selection_windows$chr[i]
      start = Balancing_selection_windows$Start[i]
      end = Balancing_selection_windows$End[i]
    }
  }
}

list_windows <- do.call(rbind, list_windows) %>% as_tibble() %>% mutate(V2 = V2-winsize, V3 = V3 +winsize)

colnames(list_windows) <- c("chr","Start","End")

list_windows$Method <- NA
print(list_windows$Method)

for(win in 1:nrow(list_windows)){
        RDA_contrib <- Res %>% filter(RDA_P_val <= 0.05) %>%
                filter(chr == list_windows$chr[win],
                       Start >=  list_windows$Start[win],
                       End <=  list_windows$End[win]) %>%
                nrow()

        Fst_contrib <- Res %>% filter(Pval_ratio <= 0.05) %>%
                filter(chr == list_windows$chr[win],
                       Start >=  list_windows$Start[win],
                       End <=  list_windows$End[win]) %>%
                nrow()


        if(RDA_contrib & Fst_contrib){
                list_windows$Method[win] <- "Shared"}
        else if(RDA_contrib){
                list_windows$Method[win] <- "RDA"}
        else if(Fst_contrib){
                list_windows$Method[win] <- "Cummulative_Fst"}

}




write.table(list_windows, "11_Res_permuts/windows_balancing_selection",
            row.names =  FALSE,
            sep ="\t",
            quote = FALSE,col.names = T)
                                                                                                    
