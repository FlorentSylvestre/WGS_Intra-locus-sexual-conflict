library(data.table)
library(tidyverse)
library(ggnewscale)


freq_F = read.table("15_exploring_SNP-by-SNP/All_two_method_female.frq", h=F, skip = 1)
freq_M = read.table("15_exploring_SNP-by-SNP/All_two_method_male.frq", h=F, skip=1)

colnames(freq_F) <- c("CHR","POS","ploidy","n","p","q")
colnames(freq_M) <- c("CHR","POS","ploidy","n","p","q")

signif <- fread("12_Signif_res/SNP_significatifs",h=T)
shared <- unlist(apply(signif[,c(7,8,9)],1,sum))
signif <- signif[which(shared >=1),] %>% mutate(color_line = shared[shared >0])

dat <- rbind(data.frame(freq_M, "sex" = "M"), data.frame(freq_F, "sex"= "F"))
dat$color_line <- 0


dat$SNP_by_SNP <-FALSE
dat$Cumulative_Fst <-FALSE
dat$RDA <-FALSE

for(i in 1:nrow(dat)){

        dat$color_line[i] <- signif$color_line[paste(signif$scaf, signif$pos, sep = "") == paste(dat$CHR[i], dat$POS[i], sep ="")]
        dat$SNP_by_SNP[i] <- signif$SNP_by_SNP[paste(signif$scaf, signif$pos, sep = "") == paste(dat$CHR[i], dat$POS[i], sep ="")]
        dat$Cumulative_Fst[i]<- signif$Cummulative_Fst[paste(signif$scaf, signif$pos, sep = "") == paste(dat$CHR[i], dat$POS[i], sep ="")]
        dat$RDA[i] <- signif$RDA[paste(signif$scaf, signif$pos, sep = "") == paste(dat$CHR[i], dat$POS[i], sep ="")]
}

dat <- dat[order(dat$color_line, decreasing =T),]
print(dat)
p <- dat %>% mutate( color_line = factor(color_line, levels = c(1,2,3))) %>%
        pivot_longer(cols = c(SNP_by_SNP,Cumulative_Fst,RDA), names_to = "Metric", values_to = "values") %>%
        filter(values == TRUE) %>%
        mutate(Metric = factor(Metric, levels = c("SNP_by_SNP", "Cumulative_Fst","RDA"))) %>%
        ggplot(., aes(x = sex, y = p, group = paste(CHR,POS, sep = ""))) +
        geom_point() +

        geom_line(aes(color = color_line)) +
        scale_color_manual(values = c("grey","black", "red")) +
        labs(color = "Shared by:")+
        theme_bw() +
        facet_grid(~Metric)
ggsave("./99_plots/Figure_S6.png",
       p,
       device = "png",
       dpi = 700,
              width = 10, height = 10)

ggsave("./99_plots/Figure_S6.pdf",
       p,
       device = "pdf",
       dpi = 700,
       width = 10, height = 10)

                                        
