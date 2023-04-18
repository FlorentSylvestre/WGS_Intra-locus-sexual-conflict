library(tidyverse)
library(data.table)
#SNP for GO Enrichment en shared between methods

Signif_SNP <- fread("12_Signif_res/SNP_significatifs", h =T)

Signif_SNP$nsnps <- rowSums(Signif_SNP[,c(7,8,9)])

Total <- Signif_SNP %>% filter(nsnps >0) %>%
        select(scaf, pos)
Shared <- Signif_SNP %>% filter(nsnps >1) %>%
        select(scaf, pos)



write.table(Total,
        "97_SNP_lists/All_signif_snps",
        row.names = F,
        quote = F,
        col.names = F,
        sep = "\t")

write.table(Shared,
        "97_SNP_lists/Shared_signif_snps",
        row.names = F,
        quote = F,
        col.names = F,
        sep = "\t")
