library(data.table)
library(tidyverse)


#importing datasets:
Res <- fread("06_Res_permuts/Res_permut.txt")
Per_SNP_Fst <- fread("08_Fst_perm/_fst_res_per_snps.txt")[,c(3:10)] %>% as.data.frame()
RDA_contrib <- fread("09_RDA/RDA_OUTPUT_RDA.contrib") %>%
        filter(Scale == "PER_WIN") %>%
        as.data.frame()

#Defining significativity thresholds:

#For windows:

Cummulative_Fst_Pval = 0.05

RDA_Pval = 0.05
Tajima <- quantile(Res$Tajima, 0.99)
###Output Significant windows:
#Adding signif infos to Res file:
Res$Signif_RDA <- Res$RDA_P_val <= RDA_Pval & Res$Tajima >= Tajima

Res$Signif_Cummulative_Fst <- Res$Pval_ratio <= Cummulative_Fst_Pval & Res$Tajima >= Tajima

##Output_significant SNP for each methods
Signif_SNP <- Per_SNP_Fst[,c(1,2,7,8)]
Signif_SNP$SNPcode <- paste(Signif_SNP[,1], Signif_SNP[,2], sep = "_")
Signif_SNP$RDA <- rep(FALSE, nrow(Signif_SNP))
Signif_SNP$Cummulative_Fst <- rep(FALSE, nrow(Signif_SNP))

#Counting number of signif SNP in each windows:
Res$RDA_count_Signif <- rep(0,nrow(Res))
Res$Cummulative_Fst_count_signif <- rep(0,nrow(Res))
FDR_Cummulative_windows <- list()
FDR_RDA <- list()

#Which snps are significant in a window:
for( win in 1:nrow(Res)){
        winchr = Res$chr[win]
        winstart = Res$Start[win]
        winend = Res$End[win]
        #RDA:

        if(Res$Signif_RDA[win]){
        winRDAcontrib = RDA_contrib %>%
                 filter(chr == winchr,
                        Start == winstart,
                        End == winend)
        signif_snps <- winRDAcontrib %>% pull(Contribution) %>%
                (function(x){lims = mean(x) + c(-1,1)*3*sd(x); return(which(x < lims[1] | x > lims[2]))})

        if(length(signif_snps) == 0){
                next()
        }
        Zcontrib <- (winRDAcontrib$Contribution - mean(winRDAcontrib$Contribution))/sd(winRDAcontrib$Contribution)
        pvalcontrib <- 2*pnorm(abs(Zcontrib),lower.tail = FALSE)
        FDR_RDA_value <- max(p.adjust(pvalcontrib, "BH")[signif_snps])
        #print(c(paste(winchr,winstart,winend, sep ="_"), FDR_RDA))
        print(FDR_RDA_value)
        print(signif_snps)
        FDR_RDA[[paste(win)]] <-  c(paste(winchr,winstart,winend, sep ="_"), FDR_RDA_value)
        Signif_SNP$RDA[Signif_SNP$SNPcode %in% paste(winRDAcontrib$chr[signif_snps], winRDAcontrib$Pos[signif_snps], sep ="_")] <- TRUE
        Res$RDA_count_Signif[win] <- length(signif_snps)


        if(Signif_SNP$RDA[1]){
                #print(winchr,winstart,winend)
                #print(winRDAcontrib[signif_snps,])
        }
        }

        #Fst
        if(Res$Signif_Cummulative_Fst[win]){
        Fst_subset <- Per_SNP_Fst %>%
                filter(scaf == winchr,
                       pos >= winstart,
                       pos <= winend)

        signif_snps <- which((1-Fst_subset$quantile_observed) <= 0.001)
        Signif_SNP$Cummulative_Fst[Signif_SNP$SNPcode %in% paste(Fst_subset$scaf[signif_snps], Fst_subset$pos[signif_snps], sep = "_")] <- TRUE

        max_signif_pval <- which.min(Fst_subset$quantile_observed[signif_snps])
        fdr <- p.adjust((1-Fst_subset$quantile_observed), "BH")[signif_snps[max_signif_pval]]
#       print(fdr)

        FDR_Cummulative_windows[[paste(win)]] <- c(paste(winchr,winstart,winend, sep ="_"), fdr)
        Res$Cummulative_Fst_count_signif[win] <- length(signif_snps)
        }
}


fwrite(Signif_SNP, "13_SNP_Signif_high_Tajima/SNP_significatifs", row.names = FALSE,  quote= FALSE, sep = "\t")
fwrite(as.data.frame(do.call(rbind,FDR_Cummulative_windows)), "13_SNP_Signif_high_Tajima/FDR_for_SNPs_in_Cummulative_Fst_windows", row.names = FALSE, quote= F, sep = "\t")

fwrite(as.data.frame(do.call(rbind,FDR_RDA)), "13_SNP_Signif_high_Tajima/FDR_for_SNPs_in_RDA", row.names = FALSE, quote= F, sep = "\t")
