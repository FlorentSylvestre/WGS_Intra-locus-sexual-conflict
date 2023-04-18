#####per windows
winsize = 250000
step = 50000

#function

roh_density_plot <- function(X1, metric,GROUP, R = "none", breaks = waiver(), limits = waiver(), guide = "none", wilcox = TRUE){
  color_real_dat <- c("black","black")  ## add "grey", first if using the subsample
  if(R != "none"){
    X1 %<>% filter(roh_group == R)
  }

  if(!("roh_group" %in% colnames(X1))){
    X1 %<>% select(Tajima,pi100,r2, Signif_RDA,Signif_Cummulative_Fst,shared) }else{
      X1 %<>% select(Tajima,pi100,r2, Signif_RDA,Signif_Cummulative_Fst,roh_group,shared) }

  type <- rep(NA, nrow(X1))
  type[(!X1$Signif_RDA & !X1$Signif_Cummulative_Fst)] <- "Non Significative"
  type[(X1$Signif_RDA | X1$Signif_Cummulative_Fst)] <- "Significative"
  Iter <- rep(NA, nrow(X1))
  Iter[(!X1$Signif_RDA & !X1$Signif_Cummulative_Fst)] <- 1002
  Iter[(X1$Signif_RDA | X1$Signif_Cummulative_Fst)] <- 1001

  X1 %<>% mutate("type" = type,  Iter = Iter) %>% filter(type == "Non Significative" | get(GROUP)) %>%
    select(-Signif_RDA,-Signif_Cummulative_Fst) %>% #,-GROUP_FST_SNP) %>%
    relocate(type,.after = Iter)
  print(X1)


 p= X1 %>%
    ggplot(aes(x = get(metric)#,
               #group = Iter
    )) +
    stat_density(aes(colour = type,lty = type ),
                           geom="line",position="identity") +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, .15)),labels = function(x) {round(x, digits = 2)})+
    scale_x_continuous(expand = c(0,0.05),breaks = breaks, limits = limits ) +
    scale_color_manual(values = color_real_dat) +
    scale_linetype_manual(values=c("solid", "dashed")) +
    labs(x ="",y ="") +
    guides(color = guide, lty = "none") +
    theme(plot.margin = margin(b=0, t=0),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          axis.line = element_line(color='black'))

  if(wilcox & length(get(metric,X1[X1$Iter ==1002, ])) & length(get(metric,X1[X1$Iter ==1001, ]))){
  max_y = ggplot_build(p)$data %>% do.call(rbind,.) %>% pull(y) %>% max()
  max_x = max(get(metric,X1))
  wilcox_pav = wilcox.test(y= get(metric,X1[X1$Iter ==1002, ]), x= get(metric,X1[X1$Iter ==1001, ]) )
  n=   length(get(metric,X1[X1$Iter ==1001, ]))
    p=p + annotate(geom = 'text', x= (max_x-0.2), y = (max_y -0.1),label = signif(wilcox_pav$`p.value`,3) ) +
      annotate(geom = 'text', x= (max_x-0.2), y = (max_y -0.2),label = paste("n =", n, sep = " "))
 

  }


    return(p)
}

library(patchwork)
library(data.table)
library(tidyr)
library(magrittr)
library(grid)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggvenn)
library(RColorBrewer)
###loading results
Res <- fread( "./12_Signif_res/Res_permut_signif.txt", header =T)

Res$pi100 <- Res$pi*100
Res$Ratio_Fst_P_value[Res$Pval_ratio ==0] <- 1/20000
#Per Snp Fst and Signif:
SNP_Signif <- fread("12_Signif_res/SNP_significatifs", h = T)

##Conv table
chrconv <- data.frame("ARABE" = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21),
                      "ROMAIN" = c("chrI","chrII","chrIII","chrIV","chrV",
                                   "chrVI", "chrVII","chrVIII", "chrIX","chrX","chrXI","chrXII",
                                   "chrXIII","chrXIV","chrXV","chrXVI","chrXVII",
                                   "chrXVIII", "chrXX", "chrXXI"))


###PLOT Figure 2:
#############MAHNATAN PLOT https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
#SNPs:
color_group_SNP <- function(SNP_signif){
        if(SNP_signif$Non_signif){
                if(SNP_signif$chr %%2){return("Non significant")}else{return("Non significant_B")}}
        if(SNP_signif$RDA){
                if(SNP_signif$SNP_by_SNP){
                        if(SNP_signif$Cummulative_Fst){return("All")}
                        else{return("RDA & SNP-by-SNP")}}
                if(SNP_signif$Cummulative_Fst){return("RDA & Cumulative Fst")}
                return("RDA")}

        if(SNP_signif$Cummulative_Fst){
                if(SNP_signif$SNP_by_SNP){return("Cumulative Fst & SNP-by-SNP")}
                return("Cumulative Fst")}

        return("SNP-by-SNP")
}

data_cum_bp <- SNP_Signif %>%
        select(scaf, pos) %>%
        group_by(scaf) %>%
        summarise(max_bp = max(pos)) %>%
        mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
        select(scaf, bp_add) %>%
        rename(chr = scaf)

SNP_Signif <- SNP_Signif %>% rename(chr = scaf) %>%
        inner_join(data_cum_bp, by = "chr") %>%
        mutate(bp_cum = pos + bp_add) %>%
        select(-bp_add)

gSNP <- rep(NA, nrow(SNP_Signif))
for(i in 1:nrow(SNP_Signif)){
        gSNP[i] <- color_group_SNP(SNP_Signif[i,])}
print(gSNP[i])


gSNP <- factor(gSNP, levels = c("Non significant", "Non significant_B",
                                "RDA & SNP-by-SNP", "RDA & Cumulative Fst",
                                "RDA", "Cumulative Fst & SNP-by-SNP",
                                "Cumulative Fst", "SNP-by-SNP", "All"))
#gSNP <- rowSums(SNP_Signif[,c(7,8,9)])

#gSNP <- factor(gSNP, levels = c("Non significant", "Non significant B",
#                               "1", "2", "3"))



axis_set_SNP <- SNP_Signif %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum))
table(gSNP)
SNP_Signif$COLORgroup <- gSNP
SNP_Signif <- SNP_Signif[order(as.numeric(gSNP)),]

SNP <- SNP_Signif %>%
  ggplot(aes(y = fst.obs, x = bp_cum, color = COLORgroup)) +
  rasterize(geom_point(size = 0.1), dpi = 300) +
  scale_x_continuous(label =chrconv$ROMAIN, breaks = axis_set_SNP$center,expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0,0), limits = range(SNP_Signif$fst.obs))  +
  scale_color_manual(values = c("grey", "lightgrey", brewer.pal(n = 8, name = 'Dark2')[-6])) +
#  scale_color_manual(values = c("grey", "lightgrey", "#56B4E9", "#009E73", "black")) +
  geom_rect(data = data.frame("XM" =  SNP_Signif %>% group_by(chr) %>% summarise(x = min(bp_cum)) %>% pull(x),
                              "XF"= SNP_Signif %>% group_by(chr) %>% summarise(x = max(bp_cum))%>% pull(x),
                              "CHROM" = c(1:21)[-19]),
            aes(xmin = XM,
                xmax =XF,
                fill =as.factor(CHROM),
                ymin = range(SNP_Signif$fst.obs)[1],
                ymax = range(SNP_Signif$fst.obs)[2]),
            alpha = 0.3,
            inherit.aes = FALSE) +
  scale_fill_manual(values = rep(c("grey","white"),10)) +
 theme_bw() +
  guides(fill = "none" ) +
  labs(y = "Intersex Fst",
       color = 'Legend :') +
  ggtitle("Significant SNPs")+
  theme(legend.position =  'bottom')


SNP_by_SNP_pval_hist <- SNP_Signif %>%
        ggplot(., aes(x= (1-quantile_observed))) +
               geom_histogram(binwidth = 0.05)+
               theme_bw()

Cumulative_Fst_Pval_hist <- Res%>%
        filter(Nsignif >0) %>%
        ggplot(.,aes(x = Pval_ratio)) +
                geom_histogram(binwidth = 0.05) +
                theme_bw()
RDA_pval_hist <- Res%>%
        ggplot(.,aes(x = RDA_P_val)) +
                geom_histogram(binwidth = 0.05) +
                theme_bw()


ggsave("99_plots//SNP_by_SNP_distrib.png", SNP_by_SNP_pval_hist, device = "png", dpi = 300)
ggsave("99_plots/Cumulative_Fst_Pval_distrib.png", Cumulative_Fst_Pval_hist, device = "png", dpi =300)
ggsave("99_plots/RDA_Pval_distrib.png", Cumulative_Fst_Pval_hist, device = "png", dpi =300)


ggsave("99_plots/SNP_by_SNP_distrib.pdf",SNP_by_SNP_pval_hist, device = "pdf", dpi =300)
ggsave("99_plots/Cumulative_Fst_Pval_distrib.pdf", Cumulative_Fst_Pval_hist, device = "pdf", dpi =300)
ggsave("99_plots/RDA_Pval_distrib.pdf", RDA_pval_hist, device = "pdf", dpi =300)

##Shared windows
Res$shared <- (Res$Signif_RDA & Res$Signif_Cummulative_Fst)

#windows
color_group_WIN =function(signif,chr,shared){

  if(signif){
    if(shared){return("shared")}
    return("Significant")}
  else{
    if( chr %%2== 0){
      return("pair")}else{return("odd")}
  }
}



data_cum_win <- Res %>% select(chr, Start, End) %>%
  mutate("POS" = (Start + End)/2) %>%
  group_by(chr) %>%
  summarise(max_bp = max(POS)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(chr, bp_add)

Res <- Res %>%
  inner_join(data_cum_win, by = "chr") %>%
  mutate(bp_cum = ((Start + End)/2) + bp_add) %>%
  select(-bp_add)

gFst = rep(NA, nrow(Res))
gRDA = rep(NA, nrow(Res))
for(i in 1:nrow(Res)){
  gFst[i] <- color_group_WIN(Res$Signif_Cummulative_Fst[i],Res$chr[i], Res$shared[i])
  gRDA[i] <- color_group_WIN(Res$Signif_RDA[i],Res$chr[i],Res$shared[i])
}


#Axis parameters
axis_set_WIN <- Res %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum))

#Signif Line for windows
Threshold_RDA <- -log(0.01)
Threshold_FST <- -log(0.01)

##Shared windows
range(Res$bp_cum)
print(range(-log(Res$RDA_P_val)))
print("Manhatan RDA")
RDA <- Res %>%
  ggplot(aes(y = -log(RDA_P_val), x = bp_cum, color =gRDA)) +
  geom_point(alpha = as.factor(gRDA), size = 0.8)  +
  scale_x_continuous(label =chrconv$ROMAIN, breaks = axis_set_WIN$center,expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0,0), limits = range(-log(Res$RDA_P_val))) +
  scale_color_manual(values = c("grey", "lightgrey", "red", "blue"),
                     labels = c('Non Significant',
                                'Non Significant',
                                'Significant windows',"shared")) +
 # #scale_alpha_manual(values = c(1, 0.75, 0.75, 1)) +
  geom_rect(data = data.frame("XM" =  Res %>% group_by(chr) %>% summarise(x = min(bp_cum)) %>% pull(x),
                              "XF"= Res %>% group_by(chr) %>% summarise(x = max(bp_cum))%>% pull(x),
                              "CHROM" = c(1:21)[-19]),
            aes(xmin = XM,
                xmax =XF,
                fill =as.factor(CHROM),
                ymin = range(-log(Res$RDA_P_val))[1],
                ymax = range(-log(Res$RDA_P_val))[2]),
            alpha = 0.3,
            inherit.aes = FALSE) +
  scale_fill_manual(values = rep(c("grey","white"),10)) +
  theme_bw() +
  guides(fill = "none" ) +
  labs(x = "none",
      y = "-log(p-value)",
      color = 'Legend :') +
  geom_hline(yintercept = -log(0.01)) +
  ggtitle("RDA significative windows")+
  annotate(geom = 'text', x = 8000000, y = -log(0.01) - 0.3, label = 'p-value 0.01', size = 3 ) +
  theme(legend.position =  'bottom',
        axis.text.x = element_blank(),
        axis.title.x=element_blank()
)


print("Manhatan Cummulative")


Cummulative_FST <- Res %>%
  ggplot(aes(y = -log(Pval_ratio), x = bp_cum, color =gFst)) +
  geom_point(alpha = as.factor(gFst), size = 0.8) +
  scale_x_continuous(label =chrconv$ROMAIN, breaks = axis_set_WIN$center,expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0,0), limits = range(-log(Res$Ratio_Fst_P_value))) +
  scale_color_manual(values = c( "grey", "lightgrey", "red", "blue"),
                     labels = c('Non Significant',
                                'Non Significant', "shared",
                                'Significant windows')) +
  scale_alpha_manual(values = c(1, 0.75, 0.75, 1)) +
  geom_rect(data = data.frame("XM" =  Res %>% group_by(chr) %>% summarise(x = min(bp_cum)) %>% pull(x),
                              "XF"= Res %>% group_by(chr) %>% summarise(x = max(bp_cum))%>% pull(x),
                              "CHROM" = c(1:21)[-19]),
            aes(xmin = XM,
                xmax =XF,
                fill =as.factor(CHROM),
                ymin = range(-log(Res$Pval_ratio))[1],
                ymax = range(-log(Res$Pval_ratio))[2]),
            alpha = 0.3,
            inherit.aes = FALSE) +
  scale_fill_manual(values = rep(c("grey","white"),10)) +
  theme_bw() +
  theme(legend.position =  'bottom',
        axis.text.x = element_blank(),
        axis.title.x=element_blank()) +
  guides(fill = "none", color = "none") +
  labs(x = "none",
       y = "-log(p-value)",
       color = 'Legend :') +
  geom_hline(yintercept = -log(0.01)) +
  ggtitle("Cumulative Fst significative windows")+
  annotate(geom = 'text', x = 8000000, y = -log(0.01) - 0.005, label = 'p-value = 0.01', size = 3 )


print("Global Manhattan")
p <- ( Cummulative_FST /RDA /SNP) +
theme(legend.position = 'bottom', plot.title  = element_text(size = 10 ))

print("Saving plot")


ggsave("./99_plots/Manhattan_plot.png",
       p,
       device = "png",
       dpi = 300,
       width = 16.9, height = 11.38, unit = 'cm')
ggsave("./99_plots/Manhattan_plot.pdf",
       p,
       device = "pdf",
       dpi = 300,
      width = 16.9, height = 11.38, unit = 'cm')


###Figure 3: Venn diagrams and Fst distrib per SNP per methods
    #####Venn diagrams of overlap between methods
dat_ven_win <- list("RDA" = paste(Res$chr[Res$Signif_RDA], Res$Start[Res$Signif_RDA], sep = "_"),
                    "Cummulative Fst" = paste(Res$chr[Res$Signif_Cummulative_Fst], Res$Start[Res$Signif_Cummulative_Fst], sep = "_"))#,
venn_win <- ggvenn(dat_ven_win, fill_color = c("white","white"))


ggsave("./99_plots/Venn_windows.png",
       plot = venn_win,
       device = "png",
       dpi = 700,
       width = 10, height = 10)
ggsave("./99_plots/Venn_windows.pdf",
       plot = venn_win,
       device = "pdf",
       dpi = 700,
       width = 10, height = 10)

#Venn diagrams of overlp between SNPs:

dat_ven_snp <- list("RDA" = paste(SNP_Signif$chr[SNP_Signif$RDA], SNP_Signif$pos[SNP_Signif$RDA], sep = "_"),
                    "Cummulative_Fst" = paste(SNP_Signif$chr[SNP_Signif$Cummulative_Fst], SNP_Signif$pos[SNP_Signif$Cummulative_Fst], sep = "_"),
                    "SNP_by_SNP" = paste(SNP_Signif$chr[SNP_Signif$SNP_by_SNP], SNP_Signif$pos[SNP_Signif$SNP_by_SNP], sep = "_"))

venn_snp <- ggvenn(dat_ven_snp, fill_color = c("white","white","white"))

ggsave("./99_plots/Venn_snp.png",
       plot = venn_snp,
       device = "png",
       dpi = 700,
       width = 10, height = 10)
ggsave("./99_plots/Venn_snp.pdf",
       plot = venn_snp,
       device = "pdf",
       dpi = 700,
       width = 10, height = 10)
 
 ############################## FST PER SNP Effect size distribution:
###wilcoxon test:
#NS - RDA
print(head(SNP_Signif))
NS_fst <- SNP_Signif %>% filter( Non_signif) %>% pull(fst.obs)
Rda_fst <- SNP_Signif %>% filter(RDA) %>% pull(fst.obs)
N_RDA_wilcox <- wilcox.test(NS_fst, y = Rda_fst)

#NS- FSTWIN
Cumulative_Fst <- SNP_Signif %>% filter(Cummulative_Fst) %>% pull(fst.obs)
NS_Cumulative_Fst_wilcox <- wilcox.test(NS_fst, y = Cumulative_Fst)

#NS- FST SNP
SNP_by_SNP_fst <- SNP_Signif %>% filter(SNP_by_SNP) %>% pull(fst.obs)
NS_SNP_by_SNP_wilcox <- wilcox.test(NS_fst, y = SNP_by_SNP_fst)

#RDA- FST WIN
RDA_FST_wilcox <- wilcox.test(Rda_fst, y = Cumulative_Fst)

#RDA SNP_by_SNP
RDA_SNP_by_SNP_wilcox <- wilcox.test(Rda_fst, y = SNP_by_SNP_fst)

#Cummulative_Fst SNP_by_SNP
SNP_by_SNP_Cumulative_FST_wilcox <- wilcox.test(Cumulative_Fst, y = SNP_by_SNP_fst)


print("NON-RDA")
print(N_RDA_wilcox)

print("Non- Cumulative")
print(NS_Cumulative_Fst_wilcox)

print("NON-SNP_BY_SNP")
print(NS_SNP_by_SNP_wilcox)

print("RDA CUMULATIVE")
print(RDA_FST_wilcox)

print("RDA SNP_by_SNP")
print(RDA_SNP_by_SNP_wilcox)

print("Cumulative SNP_by_SNP")
print(SNP_by_SNP_Cumulative_FST_wilcox)

##plots
BY_SNP_Fst <- SNP_Signif %>% pivot_longer(cols = c("Non_signif", "RDA", "Cummulative_Fst","SNP_by_SNP"), names_to= "Metric", values_to ="Signif") %>%
        filter(Signif) %>%
        mutate(Metric = factor(Metric, levels = c("Non_signif","RDA","Cummulative_Fst","SNP_by_SNP"))) %>%
        ggplot(.,aes(y = fst.obs, x = Metric)) +
        geom_point_rast(color = "grey") +
        geom_jitter(color = "grey") +
        geom_boxplot(alpha = 0)+
        theme_bw()


ggsave("./99_plots/BY_SNP_Fst.png",
       plot = BY_SNP_Fst,
       device = "png",
       dpi = 300,
       width = 10, height = 10)

ggsave("./99_plots/BY_SNP_Fst.pdf",
        plot = BY_SNP_Fst,
        device = "pdf",
        dpi = 300,
        width = 5, height = 5,
        unit = "cm",
        useDingbats = T)

ggsave("./99_plots/BY_SNP_Fst.tiff",
        plot = BY_SNP_Fst,
        device = "tiff",
        dpi = 300,
        width = 10, height = 10,
        unit = "cm")

ggsave("./99_plots/BY_SNP_Fst_lzw.tiff",
        plot = BY_SNP_Fst,
        device = "tiff",
        dpi = 300,
        width = 10, height = 10,
        unit = "cm",
        compression = "lzw")




###################PLOTS Diversity################################

####PLots distrib Tajima
##pval _ distrib_diversity

###PLOTS TAJ PI ROH per metric
#TAJIMA PLOTS
Tajima_limits <- c(-2.05,1,05)
Tajima_breaks <- c(-2,-1,0,1)
pi_limits <- c(0,0.02) *100
pi_breaks <- c(0,0.005,0.01,0.015) *100
roh_limits <- c(0,30)
roh_breaks <- seq(0,30,5)
print(Res$shared)
Taj_Fst_plot <- roh_density_plot(X1 = Res,metric = "Tajima",GROUP = "Signif_Cummulative_Fst", breaks = Tajima_breaks, limits = Tajima_limits )
pi_Fst_plot <- roh_density_plot( X1 = Res,metric = "pi100",GROUP = "Signif_Cummulative_Fst", breaks = pi_breaks, limits = pi_limits )
roh_Fst_plot <- roh_density_plot(X1 = Res,metric = "r2",GROUP = "Signif_Cummulative_Fst", breaks = roh_breaks, limits = roh_limits )
Taj_RDA_plot <- roh_density_plot(X1 = Res,metric = "Tajima",GROUP = "Signif_RDA", breaks = Tajima_breaks, limits = Tajima_limits )
pi_RDA_plot <- roh_density_plot( X1 = Res,metric = "pi100",GROUP = "Signif_RDA", breaks = pi_breaks, limits = pi_limits )
roh_RDA_plot <- roh_density_plot(X1 = Res,metric = "r2",GROUP = "Signif_RDA", breaks = roh_breaks, limits = roh_limits )

Taj_shared_plot <- roh_density_plot(X1 = Res,metric = "Tajima",GROUP = "shared", breaks = Tajima_breaks, limits = Tajima_limits )
pi_shared_plot <- roh_density_plot( X1 = Res,metric = "pi100",GROUP = "shared", breaks = pi_breaks, limits = pi_limits )
roh_shared_plot <- roh_density_plot(X1 = Res,metric = "r2",GROUP = "shared", breaks = roh_breaks, limits = roh_limits )


plot_diversity <- (Taj_RDA_plot + facet_grid(~"Tajima")) /Taj_Fst_plot /Taj_shared_plot |
(pi_RDA_plot + facet_grid(~"pi")) / pi_Fst_plot /pi_shared_plot|
(roh_RDA_plot + facet_grid("RDA"~ "roh"))/(roh_Fst_plot + facet_grid("Cummulative Fst"~.)) / (roh_shared_plot + facet_grid("Shared windows" ~.))



##PLOTS per ROH group TAJIMA
print("TAJ ROH")
Fst_roh_lowTAJ <- roh_density_plot(Res, "Tajima","Signif_Cummulative_Fst", "low", breaks = Tajima_breaks, limits = Tajima_limits)
Fst_roh_mediumTAJ <- roh_density_plot(Res, "Tajima","Signif_Cummulative_Fst", "medium", breaks = Tajima_breaks, limits = Tajima_limits)
Fst_roh_highTAJ <- roh_density_plot(Res, "Tajima","Signif_Cummulative_Fst", "high", breaks = Tajima_breaks, limits = Tajima_limits)

RDA_roh_lowTAJ <- roh_density_plot(Res, "Tajima","Signif_RDA", "low", breaks = Tajima_breaks, limits = Tajima_limits)
RDA_roh_mediumTAJ <- roh_density_plot(Res,metric =  "Tajima",GROUP = "Signif_RDA",R= "medium", breaks = Tajima_breaks, limits = Tajima_limits)
RDA_roh_highTAJ <- roh_density_plot(Res, "Tajima","Signif_RDA", "high", breaks = Tajima_breaks, limits = Tajima_limits)

Shared_roh_lowTAJ <- roh_density_plot(Res, "Tajima","shared", "low", breaks = Tajima_breaks, limits = Tajima_limits)
Shared_roh_mediumTAJ <- roh_density_plot(Res,metric =  "Tajima",GROUP = "shared",R= "medium", breaks = Tajima_breaks, limits = Tajima_limits)
Shared_roh_highTAJ <- roh_density_plot(Res, "Tajima","shared", "high", breaks = Tajima_breaks, limits = Tajima_limits)


Div_per_roh_TAJ <- (RDA_roh_lowTAJ + facet_grid(~"Low roh")) /Fst_roh_lowTAJ /Shared_roh_lowTAJ  |
(RDA_roh_mediumTAJ+ facet_grid(~"Medium roh"))/Fst_roh_mediumTAJ /Shared_roh_mediumTAJ |
(RDA_roh_highTAJ + facet_grid("RDA"~ "High roh"))/(Fst_roh_highTAJ + facet_grid("Cummulative Fst"~.))/(Shared_roh_highTAJ + facet_grid("Shared windows"~.))


##PLOTS per ROH group pi
Fst_roh_lowpi <- roh_density_plot(Res, "pi100","Signif_Cummulative_Fst", "low", breaks = pi_breaks, limits = pi_limits)
Fst_roh_mediumpi <- roh_density_plot(Res, "pi100","Signif_Cummulative_Fst", "medium", breaks = pi_breaks, limits = pi_limits)
Fst_roh_highpi <- roh_density_plot(Res, "pi100","Signif_Cummulative_Fst", "high", breaks = pi_breaks, limits = pi_limits)

RDA_roh_lowpi <- roh_density_plot(Res, "pi100","Signif_RDA", "low", breaks = pi_breaks, limits = pi_limits)
RDA_roh_mediumpi <- roh_density_plot(Res, "pi100","Signif_RDA", "medium", breaks = pi_breaks, limits = pi_limits)
RDA_roh_highpi <- roh_density_plot(Res, "pi100","Signif_RDA", "high", breaks = pi_breaks, limits = pi_limits)

Shared_roh_lowpi <- roh_density_plot(Res, "pi100","shared", "low", breaks = pi_breaks, limits = pi_limits)
Shared_roh_mediumpi <- roh_density_plot(Res,metric =  "pi100",GROUP = "shared",R= "medium", breaks = pi_breaks, limits = pi_limits)
Shared_roh_highpi <- roh_density_plot(Res, "pi100","shared", "high", breaks = pi_breaks, limits = pi_limits)

print("PI ROH")
Div_per_roh_pi <- (RDA_roh_lowpi + facet_grid(~"Low roh")) /Fst_roh_lowpi /Shared_roh_lowpi|
(RDA_roh_mediumpi+ facet_grid(~"Medium roh"))/Fst_roh_mediumpi /Shared_roh_mediumpi|
(RDA_roh_highpi + facet_grid("RDA"~ "High roh"))/(Fst_roh_highpi + facet_grid("Cummulative Fst"~.))/(Shared_roh_highpi + facet_grid("Shared windows" ~.))

print("A")
ggsave("./99_plots/Association_diversité_Metric.png",
       plot = plot_diversity,
       device = "png",
       dpi = 300,
       width = 17.9, height = 11.38,units = 'cm')
ggsave("./99_plots/Association_diversité_Metric.pdf",
       plot = plot_diversity,
       device = "pdf",
       dpi = 300,
       width = 17.9, height = 11.38,units = 'cm')

ggsave("./99_plots/Association_diversité_per_roh_pi.png",
       plot = Div_per_roh_pi,
       device = "png",
       dpi = 1200,
       width = 16.9, height = 11.38,units = 'cm')
ggsave("./99_plots/Association_diversité_per_roh_pi.pdf",
       plot = Div_per_roh_pi,
       device = "pdf",
      dpi = 1200,
       width = 16.9, height = 11.38,units = 'cm')

ggsave("./99_plots/Association_diversité_per_roh_Tajima.png",
       plot = Div_per_roh_TAJ,
       device = "png",
       dpi = 300,
       width = 17.9, height = 11.38,units = 'cm')
ggsave("./99_plots/Association_diversité_per_roh_Tajima.pdf",
       plot = Div_per_roh_TAJ,
       device = "pdf",
       dpi = 300,
       width = 17.9, height = 11.38,units = 'cm')


