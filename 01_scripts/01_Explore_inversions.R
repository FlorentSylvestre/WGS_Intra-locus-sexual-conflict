###########EXPLORE PUTATIVE GENOMIC REARRENGMENT USING LOCAL PCA.###############
#########To run manually in Rstudio would be moreefficient as this is really exploratory####

library(lostruct)
library(gridExtra)
library(patchwork)
library(tidyr)
library(dplyr)
library(magrittr)
library(ggsci)
library(windowscanr)
# Load gridExtra package
####functions

impute <- function(x){
  replace(x,
          which(x == -1),
          as.numeric(names(which.max(table(x[x!=-1]))))
  )}

plot_Twindows <- function(pca.matrix, Npca, Nind, x = 1, y = 2){

  dat_to_plot <- data.frame(
    "X" = unlist(pca.matrix[Npca, (7 + Nind * (x - 1)) : (Nind * x + 6)]),
    "Y" = unlist(pca.matrix[Npca, (7 + Nind * (y - 1)) : (Nind * y + 6)]),
    "var1" = round(pca.matrix[Npca, 5]/pca.matrix[Npca, 4], 2) * 100,
    "var2" = round(pca_matrix[Npca, 6]/pca.matrix[Npca, 4], 2) * 100, # % of variance explained by PC2
    "midpos_Npca" = (pca.matrix[Npca, 2] + pca.matrix[Npca, 3]) / 2)

  window_Npca<-paste("Windows n°:", Npca, "Chromosome :", pca_matrix[Npca, 1], "Position :", dat_to_plot$midpos_Npca , sep=" ") #paste the name of CHR and the midposition

  p<- ggplot(dat_to_plot, aes(x = X, y = Y)) +
    labs(x = paste("PC1", dat_to_plot$var1, "%"),
           y = paste("PC2", dat_to_plot$var2, "%"),
           title = window_Npca
           ) +
    geom_point() + theme_bw()

  return(p)
}
##plot mds. Add zoom in specific chrom?
mdsplot <- function(mds.matrix, mds_list = c("mds1","mds2"), manhattan = F, chr = NULL, legend = "none"){

  if( manhattan == F & length(mds_list) ==2){

        p <- ggplot(mds.matrix,
               aes(x = eval(parse(text = mds_list[1])) ,
                   y = eval(parse(text = mds_list[2])),
                   colour=chrom)) +
              geom_point() +
            theme_bw() +
              labs( x = mds_list[1], y = mds_list[2])+
              theme(legend.position = legend,
                    legend.key.size = unit(0.1,"cm"),legend.key.width = unit(0.5,"cm"))


  }else if(manhattan == T & length(mds_list) == 1){
    if(is.null(chr)){
      p <- ggplot(mds.matrix,
           aes(x=midpos,
               y = eval(parse(text = mds_list[1]))/10^6,
               colour=chrom))+
      geom_point()+
      facet_grid(cols = vars(chrom), scales = "free_x", space="free_x", switch = "x")+
      theme_bw()+
      theme(  axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_line(size = 0),
              legend.position = legend) +

      labs(x = "CHROM", y = mds_list[1])
    }else{
      sub.mds.matrix <- mds.matrix[mds.matrix$chrom == chr,]

      p <- ggplot(sub.mds.matrix,
             aes(x=midpos/10^6,
                 y = eval(parse(text = mds_list[1])))) +
        geom_point(color="blue")  +
        theme_bw() +
        labs(x = "Windows center (Mb)",
             y = mds_list[1]) +
        ggtitle(paste("chromosome ", chr,":", sep =""))+
        theme(legend.position = legend)

    }
  }
  return(p)
}

##setup variable:
color_pal <- c("#a6cee3","#1f78b4", "#b2df8a",  "#33a02c")

############################################################################################################################
#####################Global Exploration#####################################################################################
############################################################################################################################

snps <- vcf_windower("04_datasets/Full_genome/COV_MISS_QUAL_MAF010_filtered_99_autosomes_biallelic_SOR_sorted_renamed_NoUn.bcf",
                     size=1000,
                     type='snp',
                     sites= vcf_positions("04_datasets/Full_genome/COV_MISS_QUAL_MAF010_filtered_99_autosomes_biallelic_SOR_sorted_renamed_NoUn.bcf"))
pcs <- eigen_windows(snps,k=2)
window_pos<-region(snps)()
##reordering chr
tmp <- as.numeric(as.character(window_pos$chrom))
tmp[is.na(tmp)] <- 22

pcs <- pcs[order(tmp),]
window_pos <- window_pos[order(tmp),]
rm(tmp)
#Running mds

pca_matrix <- cbind(window_pos, pcs)
pcdist <- pc_dist(pcs,npc=2)
mds_axe <- cmdscale(pcdist, k=40)
mds_matrix <- cbind(window_pos, mds_axe)
colnames(mds_matrix) <- c("chrom", "start", "end", paste("mds", seq(1, 40), sep = ""))
mds_matrix$midpos <- (mds_matrix$start+mds_matrix$end)/2

mds_matrix$chrom <- factor(mds_matrix$chrom,
                           levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,"Un"))

###### Plotting all MDS, usefull for exploration
mds_biplot<- apply(data.frame(paste("mds", seq(1, 40, 2), sep=""),paste("mds", seq(2,40,2), sep="")),
                   1,
                   mdsplot, mds.matrix = mds_matrix,legend = "bottom")

##Just to look at some plots:
mds_biplot[[1]]
mds_biplot[[2]]
mds_biplot[[3]]


###plotting potentially interesting mds

mds_manhattan <- lapply(paste("mds", c( 1,2,3), sep = ""),
                        mdsplot, mds.matrix =mds_matrix, manhattan = T,legend = "none")

mds_manhattan


##########Focus on Regions of interest####
#####################CHR9#################

mds_biplot[[1]]
mds_manhattan[[1]]
md1_chr9 <- mdsplot(mds_matrix, "mds1", manhattan = T,chr = "9" )
md1_chr9
pca_to_explore <- lapply(c(810, which(mds_matrix[,"mds1"] > 0.5),840),
                           plot_Twindows,
                           pca.matrix = pca_matrix,
                           Nind = 99)
# Exploring windows around the chr9 interest region: start chr9 : 789
before <- plot_Twindows(pca_matrix,810,99)
middle <- plot_Twindows(pca_matrix,820,99)
after <- plot_Twindows(pca_matrix,840,99)

## delimit inversion on IX

md1_chr9 <- md1_chr9 + geom_vline( xintercept  = 3850000/10^6) + geom_vline(xintercept = 8600000/10^6)


##To run in a terminal to extract vcf specific to that region:
# bcftools view 04_datasets/Full_genome/COV_MISS_QUAL_MAF010_filtered_99_autosomes_biallelic_SOR_sorted_renamed_NoUn.vcf.gz -r 9:3850000-8600000 -Ob -o 04_datasets/inversions/INV_chr9.bcf
# vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 04_datasets/inversions/INV_chr9 --012

##Genotyping and ld in IV_chr9
popmap <- read.table("02_maps/sex/sexmap")
inv_9 <- read.table("/04_datasets/inversions/INV_chr9.012")[,-1]
inv_9_ind <- read.table("04_datasets/inversions/INV_chr9.012.indv")
inv_9_ind <- merge(inv_9_ind, popmap, sort =F)

rownames(inv_9) <- inv_9_ind[,1]
inv_9 <- apply(inv_9,2, impute)
inv_9_pca <- prcomp(inv_9, scale = F)

inertia <- inv_9_pca$sdev/sum(inv_9_pca$sdev) *100
geno_kmean<-kmeans (inv_9_pca$x[,1], c(min(inv_9_pca$x[,1]), (min(inv_9_pca$x[,1])+max(inv_9_pca$x[,1]))/2, max(inv_9_pca$x[,1]) ))

geno_kmean$betweenss/geno_kmean$totss * 100
PCA_chr9 <- ggplot(as.data.frame(inv_9_pca$x),aes(x = PC1, y = PC2,col = as.factor(geno_kmean$cluster), shape = inv_9_ind$V2)) +
  geom_point(size = 2.5) +
  theme_bw() +
  labs(color = "Group",
       shape = "Sex",
       x= paste("PC1 ", trunc(inertia[1]*100)/100, "%", sep = ""),
       y= paste("PC2 ", trunc(inertia[2]*100)/100, "%", sep="")) + guides(col = "none") +
  scale_color_manual(values = color_pal)

cluster_9 <- left_join(data.frame("V1" = names(geno_kmean$cluster), Group = geno_kmean$cluster), inv_9_ind)

write.table(cluster_9,
            "02_maps/maps_inversion/INV_9_cluster",
            row.names = F,
            col.names = F,
            sep = "\t",
            quote = F)
#Run outside to general sample list for each cluster / sex
#grep -P "\t1" 02_maps/maps_inversion/INV_9_cluster > 02_maps/maps_inversion/INV_9_grp1
#grep -P "\t2" 02_maps/maps_inversion/INV_9_cluster > 02_maps/maps_inversion/INV_9_grp2
#grep -P "\t3" 02_maps/maps_inversion/INV_9_cluster > 02_maps/maps_inversion/INV_9_grp3
#grep -P "\t1" 02_maps/maps_inversion/INV_9_cluster | grep M | cut -f1 > 02_maps/maps_inversion/INV_9_grp1_M
#grep -P "\t1" 02_maps/maps_inversion/INV_9_cluster | grep F | cut -f1 > 02_maps/maps_inversion/INV_9_grp1_F
#grep -P "\t2" 02_maps/maps_inversion/INV_9_cluster | grep M | cut -f1 > 02_maps/maps_inversion/INV_9_grp2_M
#grep -P "\t2" 02_maps/maps_inversion/INV_9_cluster | grep F | cut -f1 > 02_maps/maps_inversion/INV_9_grp2_F
#grep -P "\t3" 02_maps/maps_inversion/INV_9_cluster | grep M | cut -f1 > 02_maps/maps_inversion/INV_9_grp3_M
#grep -P "\t3" 02_maps/maps_inversion/INV_9_cluster | grep F | cut -f1 > 02_maps/maps_inversion/INV_9_grp3_F

##hardy_weinberg test

N <- table(cluster_9$Group)
p <- (N[1] + N[2]/2)/sum(N)
p_exp <- c(p^2, 2*p*(1-p), (1-p)^2)
Chisqtest_9 <- chisq.test(N, p = p_exp)

N_s <- table(cluster_9$Group, cluster_9$V2)
Chisqtest_9_sex <- chisq.test(N_s)

##Use vcftools to estimate per group and intersex metric for that inversion
#mkdir 98_metrics/inversion/chr9
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9 --keep 02_maps/maps_inversion/INV_9_grp1 --hardy
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9 --keep 02_maps/maps_inversion/INV_9_grp2 --hardy
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9 --keep 02_maps/maps_inversion/INV_9_grp3 --hardy
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9 --keep 02_maps/maps_inversion/INV_9_grp1 --site-mean-depth
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9 --keep 02_maps/maps_inversion/INV_9_grp2 --site-mean-depth
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9 --keep 02_maps/maps_inversion/INV_9_grp3 --site-mean-depth
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9_grp3_intersex --weir-fst-pop 02_maps/maps_inversion/INV_9_grp3_M --weir-fst-pop 02_maps/maps_inversion/INV_9_grp3_F --fst-window-size 10000 --fst-window-step 5000
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9_grp2_intersex --weir-fst-pop 02_maps/maps_inversion/INV_9_grp2_M --weir-fst-pop 02_maps/maps_inversion/INV_9_grp2_F --fst-window-size 10000 --fst-window-step 5000
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9_grp1_intersex --weir-fst-pop 02_maps/maps_inversion/INV_9_grp1_M --weir-fst-pop 02_maps/maps_inversion/INV_9_grp1_F --fst-window-size 10000 --fst-window-step 5000
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9_intersex --weir-fst-pop 02_maps/sex/M --weir-fst-pop 02_maps/sex/F --fst-window-size 10000 --fst-window-step 5000


###depth infile#
depth_9_1 <- fread("98_metrics/inversion/chr9/INV_9_grp1.ldepth.mean")
depth_9_2 <- fread("98_metrics/inversion/chr9/INV_9_grp2.ldepth.mean")
depth_9_3 <- fread("98_metrics/inversion/chr9/INV_9_grp3.ldepth.mean")

colnames(depth_9_1)[3] <- "meanCov_1"
colnames(depth_9_2)[3] <- "meanCov_2"
colnames(depth_9_3)[3] <- "meanCov_3"
depth_9 <- left_join(depth_9_1,depth_9_2, by = c("CHROM","POS")) %>%
  left_join(.,depth_9_3, by = c("CHROM","POS")) %>%
  select(-c(VAR_DEPTH.x,VAR_DEPTH.y,VAR_DEPTH)) %>%
  pivot_longer(!c(CHROM,POS)) %>%
  separate(.,name, into = c("tmp", "Group"), remove = TRUE, convert = F, sep = "_") %>%
  select(-tmp)
rm(depth_9_1,depth_9_2,depth_9_3)

##Het infile
het_9_1 <- fread("98_metricsinversion/chr9/INV_9_grp1.hwe")
het_9_2 <- fread("98_metrics/inversion/chr9/INV_9_grp2.hwe")
het_9_3 <- fread("98_metrics/inversion/chr9/INV_9_grp3.hwe")

colnames(het_9_1)[3] <- "COUNT_1"
colnames(het_9_2)[3] <- "COUNT_2"
colnames(het_9_3)[3] <- "COUNT_3"
het_9 <- left_join(het_9_1,het_9_2, by = c("CHR","POS")) %>%
  left_join(.,het_9_3, by = c("CHR","POS")) %>%
  select(c(CHR,POS, COUNT_1,COUNT_2,COUNT_3)) %>%
  pivot_longer(!c(CHR,POS)) %>%
  separate(.,name, into = c("tmp", "Group"), remove = TRUE, convert = F, sep = "_") %>%
  select(-tmp) %>%
  separate(value, into = c("HOM1", "HET", "HOM2"), remove = TRUE, convert = TRUE)


rm(het_9_1,het_9_2,het_9_3)

####fst infile
Fst_9_1 <- fread("98_metrics/inversion/chr9/INV_chr9_grp1_intersex.windowed.weir.fst")
Fst_9_2 <- fread("98_metrics/inversion/chr9/INV_chr9_grp2_intersex.windowed.weir.fst")
Fst_9_3 <- fread("98_metrics/inversion/chr9/INV_chr9_grp3_intersex.windowed.weir.fst")
Fst_9_global <- fread("98_metrics/inversion/chr9/INV_chr9_intersex.windowed.weir.fst")

colnames(Fst_9_1)[5] <- "Fst_1"
colnames(Fst_9_2)[5] <- "Fst_2"
colnames(Fst_9_3)[5] <- "Fst_3"
colnames(Fst_9_global)[5] <- "Fst_global"

Fst_win_9 <- left_join(Fst_9_1, Fst_9_2, by = c("CHROM","BIN_START","BIN_END")) %>%
  left_join(Fst_9_3, by = c("CHROM","BIN_START","BIN_END")) %>%
  left_join(Fst_9_global, by = c("CHROM","BIN_START","BIN_END")) %>%
  select(CHROM, BIN_START, BIN_END, Fst_1, Fst_2, Fst_3,Fst_global) %>%
  mutate("win_center" = (BIN_START + BIN_END)/2) %>%
  pivot_longer(c(Fst_1,Fst_2,Fst_3, Fst_global), names_to = "Group", values_to = "FST") %>%
  separate(Group, into = c("tmp","Group"), remove =T ) %>%
  select(-tmp)



##Windows

het_win <- het_9 %>%
  winScan(x = .,
          groups = "Group",
          position = "POS",
          value = c("HOM1","HET", "HOM2"),
          win_size = 10000,
          win_step = 5000,
          funs = "mean")

depth_win <- depth_9 %>%
  winScan(x = .,
          groups = "Group",
          position = "POS",
          value = "value",
          win_size = 10000,
          win_step = 5000,
          funs = "mean")

dat_win <- left_join(het_win,depth_win) %>%
  filter(!is.na(HOM1_mean ))
rm(het_win, depth_win)

###plots###
het_plot_9 <- dat_win %>%
  filter(!is.na(HOM1_mean)) %>%
  ggplot(aes(x = win_mid/10^6, y = HET_mean/(HOM1_mean + HET_mean + HOM2_mean), col = as.factor(Group)))+
  geom_line() +
  theme_bw() +
  labs(col = "Group", x = "Windows center (Mb)", y = "Heterozygosity") +
  guides(col = "none")

cov_plot_9 <- dat_win %>%
  filter(!is.na(HOM1_mean)) %>%
  ggplot(aes(x = win_mid/10^6, y = value_mean, col = as.factor(Group)))+
  geom_line() +
  theme_bw() +
  labs(col = "Group", x = "Windows center (Mb)", y = "Mean depth")+
  guides(col = "none")

Fst_plot_9 <- Fst_win_9 %>%
  ggplot(aes(x = win_center/10^6, y = FST, col = Group)) +
  geom_line() +
  theme_bw() +
  labs(x = "Position (Mb)", y = "Intersex Fst") + facet_grid(rows = vars(Group) ) +
  theme(  strip.background = element_blank(),
          strip.text.y = element_blank())

p_chr_9 <- (md1_chr9 + PCA_chr9)/(het_plot_9 | cov_plot_9| Fst_plot_9)+ plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect') & theme(axis.text.x = element_text(size = 9)) & scale_color_manual(values = color_pal)


####Zoom on potential peak in third group, between 6.5 and 7.5
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9_M --keep 02_maps/maps_inversion/INV_9_grp3_M --site-mean-depth
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9_F --keep 02_maps/maps_inversion/INV_9_grp3_F --site-mean-depth
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9_M --keep 02_maps/maps_inversion/INV_9_grp3_M --hardy
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 98_metrics/inversion/chr9/INV_chr9_F --keep 02_maps/maps_inversion/INV_9_grp3_F --hardy

zoom_pos <- c(6.5,7.5)

Cov_M_grp3_chr9 <- fread("98_metrics/inversion/chr9/INV_chr9_M.ldepth.mean") %>%
  select(CHROM,POS,MEAN_DEPTH) %>%
  rename(Cov_M = MEAN_DEPTH)
Cov_F_grp3_chr9 <- fread("98_metrics/inversion/chr9/INV_chr9_F.ldepth.mean") %>%
  select(CHROM,POS,MEAN_DEPTH) %>%
  rename(Cov_F = MEAN_DEPTH)

Cov_sex_grp3_chr9 <- inner_join(Cov_F_grp3_chr9, Cov_M_grp3_chr9) %>%
  pivot_longer(c(Cov_F,Cov_M), names_to = "Sex", values_to = "Cov") %>%
  winScan(x = .,
          position = "POS",
          group = "Sex",
          value = "Cov",
          win_size = 10000,
          win_step = 5000,
          funs = "mean") %>%
  filter(!is.na(Cov_mean)) %>%
  pivot_wider(values_from =Cov_mean, names_from = Sex) %>%
  mutate(Cov_Ratio = Cov_F/Cov_M)


het_9_3_F <- fread("98_metrics/inversion/chr9/INV_chr9_F.hwe")
het_9_3_M <- fread("98_metrics/inversion/chr9/INV_chr9_M.hwe")
colnames(het_9_3_F)[3] <- "COUNT_F"
colnames(het_9_3_M)[3] <- "COUNT_M"

het_9_3_sex <- left_join(het_9_3_F,het_9_3_M, by = c("CHR","POS")) %>%
  select(c(CHR,POS, COUNT_F,COUNT_M)) %>%
  pivot_longer(!c(CHR,POS)) %>%
  separate(.,name, into = c("tmp", "Group"), remove = TRUE, convert = F, sep = "_") %>%
  select(-tmp) %>%
  separate(value, into = c("HOM1", "HET", "HOM2"), remove = TRUE, convert = TRUE) %>%
  winScan(x = .,
          groups = "Group",
          position = "POS",
          value = c("HOM1","HET", "HOM2"),
          win_size = 10000,
          win_step = 5000,
          funs = "mean")%>%
  filter(!is.na(HET_mean)) %>%
  mutate(HET_mean_p = HET_mean/(HET_mean + HOM1_mean + HOM2_mean)) %>%
  select(win_mid,Group,HET_mean_p) %>%
  pivot_wider(values_from = HET_mean_p, names_from = Group) %>%
  mutate(Diff_het = F-M)




covsex9_zoom3 <- Cov_sex_grp3_chr9 %>%
  filter(win_mid/10^6 >= zoom_pos[1] & win_mid/10^6 <= zoom_pos[2] & !is.na(Cov_Ratio)) %>%
  ggplot(data = .) +
  geom_point(aes(x = win_mid/10^6,
                 y = log2(Cov_Ratio),
                          col = "group3")) +
  theme_bw() +
  labs(x = "",
       y = "Log2 of F:M mean coverage ratio")

Fst9_zoom3 <-Fst_win_9 %>%
  filter(win_center/10^6 >= zoom_pos[1] & win_center/10^6 <= zoom_pos[2] & Group == 3) %>%
  ggplot(data = .) +
  geom_point(aes(x = win_center/10^6,
                 y = FST, col = "group3") ) +
  theme_bw() +
  labs(x = "Windows center (Mb)", y = "Weir & Cockerham intersex Fst")

hetdiff9_zoom3 <- het_9_3_sex %>%
  filter(win_mid/10^6 >= zoom_pos[1] & win_mid/10^6 <= zoom_pos[2]) %>%
  ggplot(data = .) +
  geom_point(aes(x = win_mid/10^6,
                 y = Diff_het,
                 col = "group3")) +
  theme_bw() +
  labs(x = "",
       y = "F:M difference in mean heterozygosity")

zoom_chr9 <- covsex9_zoom3/hetdiff9_zoom3/Fst9_zoom3 &scale_color_manual(values = color_pal[3]) & guides(col = "none")


############################################################################################################################
#####################CHR16##################################################################################################
############################################################################################################################
 mds_biplot[[1]]
 mds_manhattan[[2]]
 md1_chr16 <- mdsplot(mds_matrix, "mds2", manhattan = T,chr = "16" )
 md1_chr16
 pca_to_explore <- lapply(c(1450, which(mds_matrix[,"mds2"] > 0.5),1465),
                          plot_Twindows,
                          pca.matrix = pca_matrix,
                          Nind = 99)
 #pca to keep

 before <- plot_Twindows(pca_matrix,1450,99)
 middle <- plot_Twindows(pca_matrix,1461,99)
 after <- plot_Twindows(pca_matrix,1465,99)

 ## delimit inversion on XVI

 md1_chr16 <- md1_chr16 + geom_vline( xintercept  = 17.1) + geom_vline(xintercept = 18.05)

# bcftools view <path to in dataset> -r 16:17100000-18050000 -Ob -o 04_datasets/inversions/INV_chr9.bcf
# vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 04_datasets/inversions/INV_chr9 --012

 ##Genotyping and ld in IV_chr9
popmap <- read.table("02_maps/sex/sexmap")
inv_16 <- read.table("04_datasets/inversions/INV_chr16.012")[,-1]
inv_16_ind <- read.table("04_datasets/inversions/INV_chr16.012.indv")
inv_16_ind <- merge(inv_9_ind, popmap, sort =F)



rownames(inv_16) <- inv_16_ind[,1]
inv_16 <- apply(inv_16,2, impute)
inv_16_pca <- prcomp(inv_16, scale = F)
inertia <- inv_16_pca$sdev/sum(inv_16_pca$sdev) *100

geno_kmean_16 <- kmeans (inv_16_pca$x[,1], c(min(inv_16_pca$x[,1]), (min(inv_16_pca$x[,1])+max(inv_16_pca$x[,1]))/2, max(inv_16_pca$x[,1]) ))

geno_kmean_16$betweenss/geno_kmean_16$totss * 100
PCA_chr16 <- ggplot(as.data.frame(inv_16_pca$x),
                    aes(x = PC1, y = PC2,
                        col = as.factor(geno_kmean_16$cluster), shape = inv_16_ind$V2)) +
   geom_point(size = 2.5) +
   theme_bw() +
   labs(color = "Kmeans group",
        shape = "Sex",
        x= paste("PC1 ", trunc(inertia[1]*100)/100, "%", sep = ""),
        y= paste("PC2 ", trunc(inertia[2]*100)/100, "%", sep="")) +
  guides(col = FALSE)+
  scale_color_manual(values = color_pal)

cluster_16 <- left_join(data.frame("V1" = names(geno_kmean_16$cluster), Group = geno_kmean_16$cluster), inv_16_ind)

 write.table(cluster_16,
             "02_maps/maps_inversion/INV_16_cluster",
             row.names = F,
             col.names = F,
             sep = "\t",
             quote = F)

#grep -P "\t1" 02_maps/maps_inversion/INV_16_cluster | cut -f1,2 > 02_maps/maps_inversion/INV_16_grp1
#grep -P "\t2" 02_maps/maps_inversion/INV_16_cluster | cut -f1,2 > 02_maps/maps_inversion/INV_16_grp2
#grep -P "\t3" 02_maps/maps_inversion/INV_16_cluster | cut -f1,2> 02_maps/maps_inversion/INV_16_grp3
#grep -P "\t1" 02_maps/maps_inversion/INV_16_cluster | grep M | cut -f1 > 02_maps/maps_inversion/INV_16_grp1_M
#grep -P "\t1" 02_maps/maps_inversion/INV_16_cluster | grep F | cut -f1 > 02_maps/maps_inversion/INV_16_grp1_F
#grep -P "\t2" 02_maps/maps_inversion/INV_16_cluster | grep M | cut -f1 > 02_maps/maps_inversion/INV_16_grp2_M
#grep -P "\t2" 02_maps/maps_inversion/INV_16_cluster | grep F | cut -f1 > 02_maps/maps_inversion/INV_16_grp2_F
#grep -P "\t3" 02_maps/maps_inversion/INV_16_cluster | grep M | cut -f1 > 02_maps/maps_inversion/INV_16_grp3_M
#grep -P "\t3" 02_maps/maps_inversion/INV_16_cluster | grep F | cut -f1 > 
/maps_inversion/INV_16_grp3_F


##hardy_weinberg test
#Global
N <- table(cluster_16$Group)
p <- (N[1] + N[2]/2)/sum(N)
p_exp <- c(p^2, 2*p*(1-p), (1-p)^2)
Chisqtest_16 <- chisq.test(N, p = p_exp)

#Sex
N_s <- table(cluster_16$Group,cluster_16$V2)
Chisqtest_16_Sex <- chisq.test(N_s)

#mkdir 98_metrics/inversion/chr16
#vcftools --bcf 04_datasets/inversions/INV_chr16.bcf --out 98_metrics/inversion/chr16/INV_chr16 --keep 02_maps/maps_inversion/INV_16_grp1 --hardy
#vcftools --bcf 04_datasets/inversions/INV_chr16.bcf --out 98_metrics/inversion/chr16/INV_chr16 --keep 02_maps/maps_inversion/INV_16_grp2 --hardy
#vcftools --bcf 04_datasets/inversions/INV_chr16.bcf --out 98_metrics/inversion/chr16/INV_chr16 --keep 02_maps/maps_inversion/INV_16_grp3 --hardy
#vcftools --bcf 04_datasets/inversions/INV_chr16.bcf --out 98_metrics/inversion/chr16/INV_chr16 --keep 02_maps/maps_inversion/INV_16_grp1 --site-mean-depth
#vcftools --bcf 04_datasets/inversions/INV_chr16.bcf --out 98_metrics/inversion/chr16/INV_chr16 --keep 02_maps/maps_inversion/INV_16_grp2 --site-mean-depth
#vcftools --bcf 04_datasets/inversions/INV_chr16.bcf --out 98_metrics/inversion/chr16/INV_chr16 --keep 02_maps/maps_inversion/INV_16_grp3 --site-mean-depth
#vcftools --bcf 04_datasets/inversions/INV_chr16.bcf --out 98_metrics/inversion/chr16/INV_chr16_grp3_intersex --weir-fst-pop 02_maps/maps_inversion/INV_16_grp3_M --weir-fst-pop 

/maps_inversion/INV_16_grp3_F --fst-window-size 10000 --fst-window-step 5000
#vcftools --bcf 04_datasets/inversions/INV_chr16.bcf --out 98_metrics/inversion/chr16/INV_chr16_grp2_intersex --weir-fst-pop 02_maps/maps_inversion/INV_16_grp2_M --weir-fst-pop 02_maps/maps_inversion/INV_16_grp2_F --fst-window-size 10000 --fst-window-step 5000
#vcftools --bcf 04_datasets/inversions/INV_chr16.bcf --out 98_metrics/inversion/chr16/INV_chr16_grp1_intersex --weir-fst-pop 02_maps/maps_inversion/INV_16_grp1_M --weir-fst-pop 02_maps/maps_inversion/INV_16_grp1_F --fst-window-size 10000 --fst-window-step 5000
#vcftools --bcf 04_datasets/inversions/INV_chr16.bcf --out 98_metrics/inversion/chr16/INV_chr16_intersex --weir-fst-pop 02_maps/sex/M --weir-fst-pop 02_maps/sex/F --fst-window-size 10000 --fst-window-step 5000


###depth infile#
depth_16_1 <- fread("98_metrics/inversion/chr16/INV_16_grp1.ldepth.mean")
depth_16_2 <- fread("98_metrics/inversion/chr16/INV_16_grp2.ldepth.mean")
depth_16_3 <- fread("98_metrics/inversion/chr16/INV_16_grp3.ldepth.mean")

colnames(depth_16_1)[3] <- "meanCov_1"
colnames(depth_16_2)[3] <- "meanCov_2"
colnames(depth_16_3)[3] <- "meanCov_3"
depth_16 <- left_join(depth_16_1,depth_16_2, by = c("CHROM","POS")) %>%
 left_join(.,depth_16_3, by = c("CHROM","POS")) %>%
 select(-c(VAR_DEPTH.x,VAR_DEPTH.y,VAR_DEPTH)) %>%
 pivot_longer(!c(CHROM,POS)) %>%
 separate(.,name, into = c("tmp", "Group"), remove = TRUE, convert = F, sep = "_") %>%
 select(-tmp)
rm(depth_16_1,depth_16_2,depth_16_3)

##Het infile
het_16_1 <- fread("98_metrics/inversion/chr16/INV_16_grp1.hwe")
het_16_2 <- fread("98_metrics/inversion/chr16/INV_16_grp2.hwe")
het_16_3 <- fread("98_metrics/inversion/chr16/INV_16_grp3.hwe")

colnames(het_16_1)[3] <- "COUNT_1"
colnames(het_16_2)[3] <- "COUNT_2"
colnames(het_16_3)[3] <- "COUNT_3"

het_16 <- left_join(het_16_1,het_16_2, by = c("CHR","POS")) %>%
   left_join(.,het_16_3, by = c("CHR","POS")) %>%
   select(c(CHR,POS, COUNT_1,COUNT_2,COUNT_3)) %>%
   pivot_longer(!c(CHR,POS)) %>%
   separate(.,name, into = c("tmp", "Group"), remove = TRUE, convert = F, sep = "_") %>%
   select(-tmp) %>%
   separate(value, into = c("HOM1", "HET", "HOM2"), remove = TRUE, convert = TRUE)


rm(het_16_1,het_16_2,het_16_3)

####fst infile
Fst_16_1 <- fread("98_metrics/inversion/chr16/INV_chr16_grp1_intersex.windowed.weir.fst")
Fst_16_2 <- fread("98_metrics/inversion/chr16/INV_chr16_grp2_intersex.windowed.weir.fst")
Fst_16_3 <- fread("98_metrics/inversion/chr16/INV_chr16_grp3_intersex.windowed.weir.fst")
Fst_16_global <- fread("98_metrics/inversion/chr16/INV_chr16_intersex.windowed.weir.fst")

colnames(Fst_16_1)[5] <- "Fst_1"
colnames(Fst_16_2)[5] <- "Fst_2"
colnames(Fst_16_3)[5] <- "Fst_3"
colnames(Fst_16_global)[5] <- "Fst_global"

Fst_win_16 <- left_join(Fst_16_1, Fst_16_2, by = c("CHROM","BIN_START","BIN_END")) %>%
 left_join(Fst_16_3, by = c("CHROM","BIN_START","BIN_END")) %>%
 left_join(Fst_16_global, by = c("CHROM","BIN_START","BIN_END")) %>%
 select(CHROM, BIN_START, BIN_END, Fst_1, Fst_2, Fst_3,Fst_global) %>%
 mutate("win_center" = (BIN_START + BIN_END)/2) %>%
 pivot_longer(c(Fst_1,Fst_2,Fst_3, Fst_global), names_to = "Group", values_to = "FST") %>%
 separate(Group, into = c("tmp","Group"), remove =T ) %>%
 select(-tmp)



##Windows

het_win_16 <- het_16 %>%
 winScan(x = .,
         groups = "Group",
         position = "POS",
         value = c("HOM1","HET", "HOM2"),
         win_size = 10000,
         win_step = 5000,
         funs = "mean")
depth_win_16 <- depth_16 %>%
   winScan(x = .,
           groups = "Group",
           position = "POS",
           value = "value",
           win_size = 10000,
           win_step = 5000,
           funs = "mean")

###plots###
het_plot_16 <- het_win_16 %>%
 filter(!is.na(HOM1_mean)) %>%
 ggplot(aes(x = win_mid/10^6, y = HET_mean/(HOM1_mean + HET_mean + HOM2_mean), col = as.factor(Group)))+
 geom_line() +
 theme_bw() +
 labs(col = "Group", x = "Windows center (Mb)", y = "Heterozygosity") +
 guides(col = "none") +
  scale_color_manual(values = color_pal)

cov_plot_16 <- depth_win_16 %>%
 filter(!is.na(value_mean)) %>%
 ggplot(aes(x = win_mid/10^6, y = value_mean, col = as.factor(Group)))+
 geom_line() +
 theme_bw() +
 labs(col = "Group", x = "Windows center (Mb)", y = "Mean depth")+
 guides(col = "none") +
    scale_color_manual(values = color_pal)

Fst_plot_16 <- Fst_win_16 %>%
 ggplot(aes(x = win_center/10^6, y = FST, col = Group)) +
 geom_line() +
 theme_bw() +
 labs(x = "Position (Mb)", y = "Intersex Fst") + facet_grid(rows = vars(Group) ) +
 theme(  strip.background = element_blank(),
         strip.text.y = element_blank()) +
    scale_color_manual(values = color_pal)

p_chr_16 <- (md1_chr16 + PCA_chr16)/(het_plot_16 | cov_plot_16| Fst_plot_16)+ plot_annotation(tag_levels = 'A') +
   plot_layout(guides = 'collect')



############################################################################################################################
#####################CHR21##################################################################################################
############################################################################################################################
mds_biplot[[2]]
mds_manhattan[[1]]
md1_chr21 <- mdsplot(mds_matrix, "mds3", manhattan = T,chr = "21" )
md1_chr21
pca_to_explore <- lapply(c(1700, which(mds_matrix[,"mds3"] < -0.5),1720),
                        plot_Twindows,
                        pca.matrix = pca_matrix,
                        Nind = 99)

md1_chr21 <- md1_chr21 + geom_vline( xintercept  =1800000/10^6) + geom_vline(xintercept = 2700000/10^6)

pca_to_explore <- lapply(c(1710, which(mds_matrix[,"mds3"] < -0.5),1720),
                        plot_Twindows,
                        pca.matrix = pca_matrix,
                        Nind = 99)
#pca to keep

before <- plot_Twindows(pca_matrix,1710,99)
middle <- plot_Twindows(pca_matrix,1713,99)
after <- plot_Twindows(pca_matrix,1720,99)

#bcftools view <path to in dataset> -r 21:1800000-2700000 -Ob -o 04_datasets/inversions/INV_chr21.bcf
#vcftools --bcf 04_datasets/inversions/INV_chr9.bcf --out 04_datasets/inversions/INV_chr9 --012

##Genotyping and ld in IV_chr9
popmap <- read.table("02_maps/sex/sexmap")
inv_21 <- read.table("04_datasets/inversions/INV_chr21.012")[,-1]
inv_21_ind <- read.table("04_datasets/inversions/INV_chr21.012.indv")
inv_21_ind <- merge(inv_21_ind, popmap, sort =F)



rownames(inv_21) <- inv_21_ind[,1]
inv_21 <- apply(inv_21,2, impute)
inv_21_pca <- prcomp(inv_21, scale = F)
inertia <- inv_21_pca$sdev/sum(inv_21_pca$sdev) *100

geno_kmean_21 <- kmeans (inv_21_pca$x[,2], c(min(inv_21_pca$x[,2]), (min(inv_21_pca$x[,2])+max(inv_21_pca$x[,2]))/2, max(inv_21_pca$x[,2]) ))
geno_kmean_21 <- kmeans (inv_21_pca$x[,1], c(min(inv_21_pca$x[,1]), max(inv_21_pca$x[,1]) ))

geno_kmean_21$betweenss/geno_kmean_21$totss * 100
PCA_chr21 <- ggplot(as.data.frame(inv_21_pca$x),aes(x = PC1, y = PC2,col = as.factor(geno_kmean_21$cluster), shape = inv_21_ind$V2)) +
 geom_point(size = 2.5) +
 theme_bw() +
 labs(color = "Group",
      shape = "Sex",
      x= paste("PC1 ", trunc(inertia[1]*100)/100, "%", sep = ""),
      y= paste("PC2 ", trunc(inertia[2]*100)/100, "%", sep="")) +
  guides(col = "none")+
  scale_color_manual(values = color_pal)

cluster_21 <- left_join(data.frame("V1" = names(geno_kmean$cluster), Group = geno_kmean$cluster), inv_21_ind)

write.table(cluster_21, "02_maps/maps_inversion/INV_21_cluster",
            row.names = F,
            col.names = F,
            sep = "\t",
            quote = F)

#grep -P "\t1" 02_maps/maps_inversion/INV_21_cluster | cut -f1,2 > 02_maps/maps_inversion/INV_21_grp1
#grep -P "\t2" 02_maps/maps_inversion/INV_21_cluster | cut -f1,2 > 02_maps/maps_inversion/INV_21_grp2
#grep -P "\t1" 02_maps/maps_inversion/INV_21_cluster | grep M | cut -f1 > 02_maps/maps_inversion/INV_21_grp1_M
#grep -P "\t1" 02_maps/maps_inversion/INV_21_cluster | grep F | cut -f1 > 02_maps/maps_inversion/INV_21_grp1_F
#grep -P "\t2" 02_maps/maps_inversion/INV_21_cluster | grep M | cut -f1 > 02_maps/maps_inversion/INV_21_grp2_M
#grep -P "\t2" 02_maps/maps_inversion/INV_21_cluster | grep F | cut -f1 > 02_maps/maps_inversion/INV_21_grp2_F


##hardy_weinberg test
#Global

#Sex
N_s <- table(cluster_21$Group,cluster_21$V2)
Chisqtest_21_Sex <- chisq.test(N_s)

#mkdir 98_metrics/inversion/chr16
#vcftools --bcf 04_datasets/inversions/INV_chr21.bcf --out 98_metrics/inversion/chr21/INV_chr21 --keep 02_maps/maps_inversion/INV_21_grp1 --hardy
#vcftools --bcf 04_datasets/inversions/INV_chr21.bcf --out 98_metrics/inversion/chr21/INV_chr21 --keep 02_maps/maps_inversion/INV_21_grp2 --hardy
#vcftools --bcf 04_datasets/inversions/INV_chr21.bcf --out 98_metrics/inversion/chr21/INV_chr21 --keep 02_maps/maps_inversion/INV_21_grp1 --site-mean-depth
#vcftools --bcf 04_datasets/inversions/INV_chr21.bcf --out 98_metrics/inversion/chr21/INV_chr21 --keep 02_maps/maps_inversion/INV_21_grp2 --site-mean-depth
#vcftools --bcf 04_datasets/inversions/INV_chr21.bcf --out 98_metrics/inversion/chr21/INV_chr21_grp2_intersex --weir-fst-pop 02_maps/maps_inversion/INV_21_grp2_M --weir-fst-pop 02_maps/maps_inversion/INV_21_grp2_F --fst-window-size 10000 --fst-window-step 5000
#vcftools --bcf 04_datasets/inversions/INV_chr21.bcf --out 98_metrics/inversion/chr21/INV_chr21_grp1_intersex --weir-fst-pop 02_maps/maps_inversion/INV_21_grp1_M --weir-fst-pop 02_maps/maps_inversion/INV_21_grp1_F --fst-window-size 10000 --fst-window-step 5000
#vcftools --bcf 04_datasets/inversions/INV_chr21.bcf --out 98_metrics/inversion/chr21/INV_chr21_intersex --weir-fst-pop 02_maps/sex/M --weir-fst-pop 02_maps/sex/F --fst-window-size 10000 --fst-window-step 5000


###depth infile#
depth_21_1 <- fread("98_metrics/inversion/chr21/INV_21_grp1.ldepth.mean")
depth_21_2 <- fread("98_metrics/inversion/chr21/INV_21_grp2.ldepth.mean")

colnames(depth_21_1)[3] <- "meanCov_1"
colnames(depth_21_2)[3] <- "meanCov_2"
depth_21 <- left_join(depth_21_1,depth_21_2, by = c("CHROM","POS")) %>%
  select(-c(VAR_DEPTH.x,VAR_DEPTH.y)) %>%
  pivot_longer(!c(CHROM,POS)) %>%
  separate(.,name, into = c("tmp", "Group"), remove = TRUE, convert = F, sep = "_") %>%
  select(-tmp)
rm(depth_21_1,depth_21_2)

##Het infile
het_21_1 <- fread("98_metrics/inversion/chr21/INV_21_grp1.hwe")
het_21_2 <- fread("98_metrics/inversion/chr21/INV_21_grp2.hwe")


colnames(het_21_1)[3] <- "COUNT_1"
colnames(het_21_2)[3] <- "COUNT_2"

het_21 <- left_join(het_21_1,het_21_2, by = c("CHR","POS"))%>%
  select(c(CHR,POS, COUNT_1,COUNT_2)) %>%
  pivot_longer(!c(CHR,POS)) %>%
  separate(.,name, into = c("tmp", "Group"), remove = TRUE, convert = F, sep = "_") %>%
  select(-tmp) %>%
  separate(value, into = c("HOM1", "HET", "HOM2"), remove = TRUE, convert = TRUE)

rm(het_21_1,het_21_2)

####fst infile
Fst_21_1 <- fread("98_metrics/inversion/chr21/INV_chr21_grp1_intersex.windowed.weir.fst")
Fst_21_2 <- fread("98_metrics/inversion/chr21/INV_chr21_grp2_intersex.windowed.weir.fst")
Fst_21_global <- fread("98_metrics/inversion/chr21/INV_chr21_intersex.windowed.weir.fst")

colnames(Fst_21_1)[5] <- "Fst_1"
colnames(Fst_21_2)[5] <- "Fst_2"
colnames(Fst_21_global)[5] <- "Fst_global"

Fst_win_21 <- left_join(Fst_21_1, Fst_21_2, by = c("CHROM","BIN_START","BIN_END")) %>%
  left_join(Fst_21_global, by = c("CHROM","BIN_START","BIN_END")) %>%
  select(CHROM, BIN_START, BIN_END, Fst_1, Fst_2,Fst_global) %>%
  mutate("win_center" = (BIN_START + BIN_END)/2) %>%
  pivot_longer(c(Fst_1,Fst_2, Fst_global), names_to = "Group", values_to = "FST") %>%
  separate(Group, into = c("tmp","Group"), remove =T ) %>%
  select(-tmp)



##Windows

het_win_21 <- het_21 %>% select(CHR,POS,Group,HET, HOM1,HOM2,) %>%
  winScan(x = .,
          groups = "Group",
          position = "POS",
          value = c("HOM1","HET", "HOM2"),
          win_size = 10000,
          win_step = 5000,
          funs = "mean")
depth_win_21 <- depth_21 %>% select(CHROM,POS,Group, value) %>%
  winScan(x = .,
          groups = "Group",
          position = "POS",
          value = "value",
          win_size = 10000,
          win_step = 5000,
          funs = "mean")

###plots###
het_plot_21 <- het_win_21 %>%
  filter(!is.na(HOM1_mean)) %>%
  ggplot(aes(x = win_mid/10^6, y = HET_mean/(HOM1_mean + HET_mean + HOM2_mean), col = as.factor(Group)))+
  geom_line() +
  theme_bw() +
  labs(col = "Group", x = "Windows center (Mb)", y = "Heterozygosity") +
  guides(col = "none") +
  scale_color_manual(values = color_pal)

cov_plot_21 <- depth_win_21 %>%
  filter(!is.na(value_mean)) %>%
  ggplot(aes(x = win_mid/10^6, y = value_mean, col = as.factor(Group)))+
  geom_line() +
  theme_bw() +
  labs(col = "Group", x = "Windows center (Mb)", y = "Mean depth")+
  guides(col = "none") +
  scale_color_manual(values = color_pal)

Fst_plot_21 <- Fst_win_21 %>%
  ggplot(aes(x = win_center/10^6, y = FST, col = Group)) +
  geom_line() +
  theme_bw() +
  labs(x = "Position (Mb)", y = "Intersex Fst") + facet_grid(rows = vars(Group) ) +
  theme(  strip.background = element_blank(),
          strip.text.y = element_blank()) +
  scale_color_manual(values = color_pal)

p_chr_21 <- (md1_chr21 + PCA_chr21)/(het_plot_21 | cov_plot_21| Fst_plot_21)+ plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect')


#Saving all plots:
ggsave(filename = "99_plots/SumPlot_INV_9.png",
       device = "png",
       plot = p_chr_9 )
ggsave(filename = "99_plots/SumPlot_INV_9.pdf",
       device = "pdf",
       plot = p_chr_9)

ggsave(filename = "99_plots/ZOOM_SumPlot_INV_9.png",
       device = "png",
       plot = p_chr_9 )
ggsave(filename = "99_plots/ZOOM_SumPlot_INV_9.pdf",
       device = "pdf",
       plot = zoom_chr9)

ggsave(filename = "99_plots/SumPlot_INV_16.png",
       device = "png",
       plot = p_chr_16 )
ggsave(filename = 99_plots/SumPlot_INV_16.pdf",
       device = "pdf",
       plot = p_chr_16 )

ggsave(filename = "99_plots/SumPlot_INV_21.png",
       device = "png",
       plot = p_chr_21 )
ggsave(filename = "99_plots/SumPlot_INV_21.pdf",
       device = "pdf",
       plot = p_chr_21 )
