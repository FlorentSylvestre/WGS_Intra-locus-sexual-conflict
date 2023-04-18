library(data.table)
library(tidyverse)
library(ggnewscale)


freq_F = read.table("15_exploring_SNP-by-SNP/Shared_two_method_female.frq", h=F, skip = 1)
freq_M = read.table("15_exploring_SNP-by-SNP/Shared_two_method_male.frq", h=F, skip=1)

colnames(freq_F) <- c("CHR","POS","ploidy","n","p","q")
colnames(freq_M) <- c("CHR","POS","ploidy","n","p","q")

signif <- fread("12_Signif_res/SNP_significatifs",h=T)
shared <- signif[apply(signif[,c(7,8,9)],1, function(x){sum(x)}) == 3, c(1,2)]
shared <- paste(shared$scaf, shared$pos, sep ="_")

dat <- rbind(data.frame(freq_M, "sex" = "M"), data.frame(freq_F, "sex"= "F"))


dat$color_line <- 2
dat$color_line[which(paste(dat$CHR,dat$POS, sep ="_") %in% shared)] <- 3
p <- dat %>% mutate(color_line = factor(color_line, levels = c(2,3))) %>%
        ggplot(., aes(x = sex, y = p, group = paste(CHR,POS, sep = ""))) +
        geom_point() +
        geom_line(aes(color = color_line)) +
        scale_color_manual(values = c("grey","black")) +
        labs(color = "Shared by:")+
        theme_bw()

ggsave("./99_plots/Figure_3A.png",
       p,
       device = "png",
       dpi = 700,
       width = 10, height = 10)

ggsave("./99_plots/Figure_3A.pdf",
       p,
       device = "pdf",
       dpi = 700,
       width = 10, height = 10)
