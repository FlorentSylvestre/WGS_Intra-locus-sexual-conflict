#usage:
#PCA_global.R <prefix> <popmap>
#prefix : prefix to  012 files; ex : test to read test.012, test.012.indv 7 test.012.pos
#popmap :: 2tab delimited colun: first = ind name, second = pop. no header
#FILE SHOULD BE IMPUTED, will not be performed in this script Use 02_impute.py to perform simple imputation

#library

library(data.table)
library(ggplot2)
library(magrittr)
if(length(commandArgs(trailingOnly = T))){
##read input file
dat <- as.data.frame(fread(paste0(commandArgs(trailingOnly = T)[1],".012"))[,-1])
ind = as.data.frame(fread(paste0(commandArgs(trailingOnly = T)[1],".012.indv"),header = F))
popmap <- read.table(commandArgs(trailingOnly = T)[2])
NPC_to_plot <- as.numeric(commandArgs(trailingOnly = T)[3])
output <- commandArgs(trailingOnly = T)[4]}else{

dat <- as.data.frame(fread("04_datasets/Full_genome/IMPUTED_NO_INV_COV_MISS_QUAL_MAF010_filtered_99_autosomes_biallelic_SOR_sorted_renamed_NoUn.012")[,-1])
ind = as.data.frame(fread("04_datasets/Full_genome/IMPUTED_NO_INV_COV_MISS_QUAL_MAF010_filtered_99_autosomes_biallelic_SOR_sorted_renamed_NoUn.012.indv",header = F))

popmap <- read.table("02_maps/sex/sexmap")
NPC_to_plot <- 12
output <- "05_PCA/Global/Global_pca"
}




##running pca
ind <- merge(ind, popmap, sort = F)
pca <- prcomp(dat,scale =F)

###plots
PCA_plot <- function(pca,grouping, pcvector = c(1,2)){
sd <- pca$sdev/sum(pca$sdev) *100
pca$x %>%
  as.data.frame() %>%
  ggplot(aes_string(x = paste0("PC",pcvector[1]),
                    y = paste0("PC",pcvector[2]))) +
  geom_point(aes(col = grouping,
                 shape = grouping))+
  stat_ellipse(aes(col =as.factor(grouping)))+

  theme_bw() +
  labs(x = paste0("PC",pcvector[1], " ", trunc(sd[pcvector[1]] *100)/100,"%"),
       y = paste0("PC",pcvector[2], " ", trunc(sd[pcvector[2]] *100)/100,"%"),
       col = "Sex")

}


List_plot <- list()
for( i in seq (1:NPC_to_plot)[c(T,F)]){
if(i != NPC_to_plot){
  List_plot[[paste(i)]] <- PCA_plot(pca,ind$V2,c(i,i+1))}else{
List_plot[[paste(i)]] <- PCA_plot(pca,ind$V2,c(i,i-1))

  }
  ggsave(paste0(output,"_",i,".png"),plot = List_plot[[paste(i)]] + scale_color_manual(values = c("red","blue")),device = "png")
  ggsave(paste0(output,"_",i,".pdf"),plot = List_plot[[paste(i)]] + scale_color_manual(values = c("red","blue")),device = "pdf")


}


#To run mannually to fine tuning of plots
PC1_2 <- List_plot[[1]] +
 scale_y_continuous(limits = c(-25, 50),
                     breaks = seq(-25,50,25)) +
  scale_x_continuous(limits = c(-50,50),
                     breaks = seq(-50,50,25)) +
  scale_color_manual(values = c("red","blue")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

PC3_4 <- List_plot[[2]] +
  scale_y_continuous(limits = c(-300, 300),
                     breaks = seq(-300,300,50)) +
  scale_x_continuous(limits = c(-75,75),
                     breaks = seq(-75,75,25)) +
   scale_color_manual(values = c("red","blue"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("99_plots/Global_PCA_1_2_Zoom.png",
       device = "png",
       plot = PC1_2,
       dpi = 300)

ggsave("99_plots/Global_PCA_1_2_Zoom.pdf",
       device = "pdf",
       plot = PC1_2,
       dpi = 300)

ggsave("99_plots/Global_PCA_3_4_Zoom.png",
       device = "png",
       plot = PC3_4,
       dpi = 300)

ggsave("99_plots/Global_PCA_3_4_Zoom.pdf",
       device = "pdf",
       plot = PC3_4,
       dpi = 300)
