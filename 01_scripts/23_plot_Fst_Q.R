library(tidyverse)
library(data.table)
library(magrittr)
library(ggnewscale)
library(patchwork)

###chr
chrconv <- data.frame("ARABE" = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21),
                      "ROMAIN" = c("chrI","chrII","chrIII","chrIV","chrV",
                                   "chrVI", "chrVII","chrVIII", "chrIX","chrX","chrXI","chrXII",
                                   "chrXIII","chrXIV","chrXV","chrXVI","chrXVII",
                                   "chrXVIII", "chrXX", "chrXXI"))
###Fst
#global
datFstpos<- data.frame(fread("08_Fst_perm/_fst_res_per_snps.txt"))
datFstglob<- t(fread("08_Fst_perm/_global_count_percentile.txt"))
colnames(datFstglob) <- seq(1,100)
colnames(datFstpos[,9:ncol(datFstpos)]) <- paste0(colnames(datFstpos[,9:ncol(datFstpos)]), sort(rep(c(1:101),2)))


#Perchr
datFstchr_perm<- fread("08_Fst_perm/perm_by_chr_ratio", header = T)
datFstchr_pval <- fread("08_Fst_perm/perm_by_chr_pval", h=T)

#function
count_Q<-function(X,Q){ apply(X,2,function(x,Q){return(sum(x>=Q,na.rm = TRUE))},Q = Q)}

#plot Globaux
#Fst_bis

Fst_999_Count = count_Q(datFstpos[,c(9:ncol(datFstpos))[c(F,T)]],0.999)
Fst_999_R = as.numeric(Fst_999_Count[1]) / as.numeric(Fst_999_Count[-1])
mean_R999 <- mean(log2(Fst_999_R))
sd_R999 <- sd(log2(Fst_999_R))

Fst_95 <- cbind(rownames(datFstglob), as.numeric(apply(datFstglob[,c(95,96,97,98,99,100)],1,sum)))
Fst_95_R <- as.numeric(Fst_95[Fst_95[,1] == "Obs",2]) / as.numeric(Fst_95[Fst_95[,1] != "Obs",2])
mean_R95 <- mean(log2(Fst_95_R))
sd_R_95 <- sd(log2(Fst_95_R))


Fst_99 <- cbind(rownames(datFstglob), datFstglob[,100])
Fst_99_R <- as.numeric(Fst_99[Fst_99[,1] == "Obs",2]) / as.numeric(Fst_99[Fst_99[,1] != "Obs",2])
mean_R99 <- mean(log2(Fst_99_R))
sd_R_99<- sd(log2(Fst_99_R))

Fst_glob <- as.data.frame(rbind(cbind("95", Fst_95_R), cbind("99", Fst_99_R),cbind("99.9",Fst_999_R)))
colnames(Fst_glob) <- c("Q", "Count")
sumst_glob <- data.frame("Mean" = c(mean_R95, mean_R99,mean_R999), "sd" = c(sd_R_95, sd_R_99,sd_R999), "Q" = c("95","99","99.9"))
Fst_glob$Count <- as.numeric(as.character(Fst_glob$Count))
  plotGlob <- Fst_glob %>%
  ggplot(aes(x = Q, y = Count)) + geom_point(col = "grey", alpha = 0.5) +
  geom_point(data = sumst_glob, mapping = aes(x = Q, y = 2^Mean), col = "red") +
  geom_segment(data = sumst_glob, mapping =  aes(x = Q, xend = Q, y = 2^(Mean + 1.96 * sd / sqrt(100)), yend =  2^(Mean - 1.96 * sd / sqrt(100))))+
  geom_hline(yintercept = 1)+
  labs(y = "Ratio of observed:permutted SNPs", x = "Fst quantiles") +
    theme_bw()


#Fst_per_chr
datFstchr <- datFstchr_perm %>%
        group_by(chr,signif_thresh) %>%
        summarise(mean = mean(log2(ratio)),
                  sd = sd(log2(ratio))) %>%
                merge(.,datFstchr_pval)

print(datFstchr$chr)
print(as.factor(datFstchr_perm$chr))
per_chr_plot <- datFstchr_perm %>%
  filter(chr != "Un") %>%
  ggplot(aes(x = as.factor(chr),
             y = ratio,
             group = as.factor(100*(1-signif_thresh)))) +
  geom_point(alpha = 0.2,
             position = position_dodge(width = 0.5),
             col = "grey") +
  geom_hline(yintercept = 1) +
  labs(y = "Ratio of observed:permutted SNPs",
       x = "Chromosomes")+
  theme_bw() +
  geom_rect(data = data.frame("XM" =  (1:20 -0.5),
                              "XF"= (1:20 +0.5),
                              "chr" = 1:20),
            aes(xmin = XM,
                xmax =XF,
                fill =as.factor(chr),
                ymin = range(datFstchr_perm$ratio)[1]-0.1 ,
                ymax = range(datFstchr_perm$ratio)[2]+0.1),
            alpha = 0.3,
            inherit.aes = FALSE)+
    scale_y_continuous(limits = range(datFstchr_perm$ratio)+c(-0.1,0.1),
                       expand = c(0, 0))+
  scale_x_discrete(expand = c(0,0))+
  scale_fill_manual(values = rep(c("grey","white"),11)[-22]) +
  theme(legend.position = c(0.075,0.7),
        legend.key=element_blank()) +
  guides(fill = "none")+
  ggnewscale::new_scale_color()+
  geom_point(data= datFstchr,
             mapping = aes(x= as.factor(chr),
                           y = 2^mean,
                           group = as.factor(100*(1-signif_thresh)),
                           col = as.factor(100*(1-signif_thresh))),
             position = position_dodge2(width = 0.5))+
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a")) +
  geom_linerange(data=datFstchr, inherit.aes = FALSE,
                 mapping = aes(x = as.factor(chr) ,
                               ymin = 2^(mean + 1.96*sd/10),
                               ymax = 2^(mean - 1.96*sd/10),
                               group = as.factor(100*(1-signif_thresh))),
                 position = position_dodge(width = 0.5)) +
                 labs(color = "Fst Quantile")
  
plot_final <- plotGlob/per_chr_plot



ggsave(filename = paste0("99_plots/Noinv_before_maf10_","Fst_permut.png"),plot = plot_final,
       device="png",
      dpi = 300,
      width = 16.9, units = 'cm')
ggsave(filename = paste0("99_plots/Noinv_before_maf10_","Fst_permut.pdf"),plot = plot_final
       , device="pdf",
       dpi = 300,
       width = 16.9, units = 'cm')


