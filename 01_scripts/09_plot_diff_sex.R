#Rscript generate the figure 1C plot


#library 
library(tidyverse)
library(data.table)
library(gridExtra)
library(patchwork)
library(ggrastr)

#Input files
dat_covhet <- fread("98_metrics/Intersex_covhet")
dat_countsnp <-  "98_metrics/intersex_snp_count_differences.txt")
dat_wilcox<- fread("98_metrics/wilcoxonRankSum_intersex_coverage")

#ordering countsnp:
dat_countsnp <- dat_countsnp[order(dat_countsnp$CHROM,dat_countsnp$POS),]


##graphic parameters
alpha = 0.5
size = 11

dat_covhet <- dat_covhet %>%
  pivot_wider(names_from = SEX, values_from = c(Het,hom, Median_cov, N))

logical = dat_wilcox$pval<=0.05
logical2 = dat_countsnp$P_fish <= 0.05

group <- rep(NA,length(logical))
for(i in 1: length(logical)){
  if(!logical[i] & !logical2[i] ){
    group[i] <- "PASS"
  }else if (logical[i] & !logical2[i]){
    group[i] <- "COV"
  }else if(!logical[i] & logical2[i]){
    group[i] <- "NSNP"
  }else{group[i] <- "BOTH"}
}
group <- factor(group, levels = c("PASS", "COV","NSNP","BOTH"))

p1 <-
  dat_covhet %>%
  ggplot(aes(y = log2(Median_cov_F/Median_cov_M), x = abs(Het_F/N_F - Het_M/N_M), col = group)) +
  rasterise(geom_point(alpha = alpha),dpi = 300) +
  scale_x_continuous(expand = c(0,0))+
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(y = "Log2 of F:M median coverage ratio",
       col = "Filter",
       x = "Heterozygosity difference between sexes") +
  theme(axis.title.y = element_text(size = size),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("#A3BBC8","#83a920","#04bbc4","#f4746a"),
                     labels = c("Pass", "Coverage Biais", "SNP count biais", "Both"))+
  geom_hline(yintercept = 0, color = 'black', linetype  = 'dashed')



ggsave(filename ="99_plots/plot_potential_duplication.png",
       device = "png",
       plot =  p1,
       dpi = 300, unit = 'mm', width = 82, height = 52)

ggsave(filename = "99_plots/plot_potential_duplication.pdf",
       device = "pdf",
       plot =  p1,
       dpi = 300, unit = 'mm', width = 82, height = 52)



fwrite(cbind(dat_covhet[group != "PASS", c("CHROM", "POS")],group[group != "PASS" ]),
       "97_SNPs_lists/Sex_suplications/BAD_Snps",
             row.names = FALSE,
             col.names = FALSE,
             sep = "\t")
fwrite(cbind(dat_covhet[group == "PASS", c("CHROM", "POS")],
             group[group == "PASS" ]),
       "97_SNPs_lists/Sex_suplications/GOOD_Snps",
       row.names = FALSE,
       col.names = FALSE,
       sep = "\t")

count_cov <- dat_covhet %>%
  filter(group %in% c("COV","BOTH")) %>%
  mutate(Ratio = Median_cov_F/Median_cov_M) %>%
  select(Ratio) %>% mutate(biais = Ratio >1) %>%
  pull(biais) %>% table
count_NSP <- dat_countsnp %>%
  filter(group %in% c("NSNP","BOTH")) %>%
  mutate(Nsnpdiff = Nf-Nm) %>%
  select(Nsnpdiff) %>% mutate(biais = Nsnpdiff > 0) %>%
  pull(biais) %>% table


sumtable <- data.frame("Females" = c(count_cov[1],count_NSP[1]),
                               "Males" = c(count_cov[2],0), "Total" = c(sum(count_cov),sum(count_NSP)))
row.names(sumtable) <- c("Coverage biais", "SNP count biais")
write.csv(sumtable,
          "97_SNPs_lists/Duplicated_summary_table.csv",
          quote = FALSE,
          row.names = TRUE)
