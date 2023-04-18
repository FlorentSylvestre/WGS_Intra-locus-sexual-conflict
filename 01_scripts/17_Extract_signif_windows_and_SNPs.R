library(data.table)
library(tidyverse)


#importing datasets:
Res <- fread("11_Res_permuts/Res_permut.txt")
Per_SNP_Fst <- fread("08_Fst_perm/fst_res_per_snps.txt")[,c(3:10)] %>% as.data.frame()
RDA_contrib <- fread("09_RDA/OUTPUT_RDA.contrib") %>%
        filter(Scale == "PER_WIN") %>%
        as.data.frame()
Fisher <- fread("07_fisher_test/fishser_test.pval", h=T)

#Defining significativity thresholds:

#For windows:

Cummulative_Fst_Pval = 0.01
Cummulative_Fst_Ratio = quantile(Res$Ratio_Fst,0.99)

RDA_Pval = 0.01
RDA_R2 = quantile(Res$RDA_R2, 0.99)

#For Snp-by-snp:

Per_SNP_Fst$qval <- p.adjust((1-Per_SNP_Fst$quantile_observed),"BH")

Threshold_SNP_by_SNP <- 0.01


###Output Significant windows:
#Adding signif infos to Res file:
print("a")
Res$Signif_RDA <- (Res$RDA_P_val <= RDA_Pval)# & Res$RDA_R2 >= RDA_R2)

Res$Signif_Cummulative_Fst <- (Res$Pval_ratio <= Cummulative_Fst_Pval)



##Output_significant SNP for each methods
Signif_SNP <- Per_SNP_Fst[,c(1,2,7,8)]
Signif_SNP$SNPcode <- paste(Signif_SNP[,1], Signif_SNP[,2], sep = "_")
Signif_SNP$Non_signif <- rep(TRUE, nrow(Signif_SNP))
Signif_SNP$RDA <- rep(FALSE, nrow(Signif_SNP))
Signif_SNP$SNP_by_SNP <- rep(FALSE, nrow(Signif_SNP))
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


        Zcontrib <- (winRDAcontrib$Contribution - mean(winRDAcontrib$Contribution))/sd(winRDAcontrib$Contribution)
        pvalcontrib <- 2*pnorm(abs(Zcontrib),lower.tail = FALSE)
        print(range(winRDAcontrib$Contribution))
        print(range(Zcontrib))
        print(range(pvalcontrib))
        if(length(signif_snps) >0){
        fdr_winrda <- max(p.adjust(pvalcontrib, "BH")[signif_snps])}
        FDR_RDA[[paste(win)]] <-  c(paste(winchr,winstart,winend, sep ="_"), fdr_winrda)

        Signif_SNP$RDA[Signif_SNP$SNPcode %in% paste(winRDAcontrib$chr[signif_snps], winRDAcontrib$Pos[signif_snps], sep ="_")] <- TRUE
        Res$RDA_count_Signif[win] <- length(signif_snps)
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

        FDR_Cummulative_windows[[paste(win)]] <- c(paste(winchr,winstart,winend, sep ="_"), fdr)
        Res$Cummulative_Fst_count_signif[win] <- length(signif_snps)
        }
}

###SNP-by-SNP significativity:
Signif_SNP$SNP_by_SNP[ p.adjust(Fisher$p_value, "BH") <= 0.05] <- TRUE


print(sum(Signif_SNP$SNP_by_SNP))
##non-Significant_SNPs:
Signif_SNP$Non_signif <- unlist(apply(Signif_SNP[,c("RDA","Cummulative_Fst", "SNP_by_SNP")],
                                      1,
                                      function(x){return(!(x[1] |x[2] |x[3]))}))



#writing results:
fwrite(Res, "12_Signif_res/Res_permut_signif.txt", row.names = FALSE,  quote= FALSE, sep = "\t")
fwrite(Signif_SNP, "12_Signif_res/SNP_significatifs", row.names = FALSE,  quote= FALSE, sep = "\t")
fwrite(as.data.frame(do.call(rbind,FDR_Cummulative_windows)), "12_Signif_res/FDR_for_SNPs_in_Cummulative_Fst_windows", row.names = FALSE, quote= F, sep = "\t")
fwrite(as.data.frame(do.call(rbind,FDR_RDA)), "12_Signif_res/FDR_for_SNPs_in_RDA", row.names = FALSE, quote= F, sep = "\t")


                                                        
