#Packages
library(patchwork)
library(data.table)
library(tidyr)
library(magrittr)
library(grid)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(ggrastr)

###loading results
Signif_SNP <- fread("13_SNP_Signif_high_Tajima/SNP_significatifs")
signif_windows <- fread("11_Res_permuts/windows_balancing_selection")
###tajima in smaler windows:

Theta_25kb <- fread("10_theta/No_max_depth.thetaswindow.pestPG") %>%
  mutate("pi" = tP/nSites) %>%
  select(Chr,WinCenter,Tajima,pi)
##Conv table
chrconv <- data.frame("ARABE" = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21),
                      "ROMAIN" = c("chrI","chrII","chrIII","chrIV","chrV",
                                   "chrVI", "chrVII","chrVIII", "chrIX","chrX","chrXI","chrXII",
                                   "chrXIII","chrXIV","chrXV","chrXVI","chrXVII",
                                   "chrXVIII", "chrXX", "chrXXI"))

print("A")
for(i in 1:nrow(Theta_25kb)){
  Theta_25kb$Chr[i]<-  chrconv$ARABE[chrconv$ROMAIN == Theta_25kb$Chr[i]]
}

#Coloring group for SNP
print("A")
set_color <- function(x){
        print(x)
        if(x[,6]){
                   if(x[,7]){ return("Both")}
                   return("RDA")}
        if(x[,7]){return("Cummulative Fst")}
        return("None")
}
print("A")
Signif_SNP$color <- rowSums(Signif_SNP[,c(6,7)])

print(Signif_SNP)

####ld/gene plots

for(metric in "Shared"){
        list_plot <- list()
        for(win in 1:nrow(signif_windows[signif_windows$Method == metric])){
                 ld <- fread(paste("15_ld/",
                            signif_windows[signif_windows$Method == metric,][win,1],
                         "_",
                        signif_windows[signif_windows$Method == metric,][win,2],
                         "-",
                        ".ldblock",sep =""))
                print(win)
                sub_theta<- Theta_25kb %>%
                filter(Chr == signif_windows[signif_windows$Method == metric,]$chr[win] &
                       WinCenter -12500 >= signif_windows[signif_windows$Method == metric,]$Start[win] &
                       WinCenter + 12500 <= signif_windows[signif_windows$Method == metric,]$End[win])
                head(Signif_SNP)
                p<- Signif_SNP %>% filter(scaf == signif_windows[signif_windows$Method == metric,]$chr[win] &
                                            pos >= signif_windows[signif_windows$Method == metric,]$Start[win] &
                                            pos <= signif_windows[signif_windows$Method == metric,]$End[win],) %>%
                                mutate(p = 1-quantile_observed) %>% pull(p)

                q <- p.adjust(p, "BH")[which(p == max(p[p<=0.001]))]
                print(q)

                p1 = plot_LD(ld[ld$Percent == 0.01,], "R2") +
                        labs( x ='Position (Mb)', y = "Position (Mb)") +
                        ggtitle(paste(chrconv$ROMAIN[chrconv$ARABE == signif_windows[signif_windows$Method == metric,]$chr[win]],
                                        ": ",
                                        format(signif_windows[signif_windows$Method == metric,]$Start[win], scientific = F),
                                        "-",
                                        format(signif_windows[ signif_windows$Method == metric,]$End[win], scientific = F),
                                        sep =""))+
                        coord_fixed() +
                        geom_vline(data = Signif_SNP %>%
                                        filter(scaf == signif_windows[signif_windows$Method == metric,]$chr[win] &
                                                pos >= signif_windows[signif_windows$Method == metric,]$Start[win] &
                                                pos <= signif_windows[signif_windows$Method == metric,]$End[win],
                                                color >0) %>%
                                        mutate(Y = max(sub_theta$Tajima)),
                                mapping = aes( xintercept=(pos/10^6), color = color), col = "red",size = 0.1, inherit.aes = FALSE)



                p3 <- ggplot(sub_theta) +
                        geom_smooth(aes(x = WinCenter/10^6, y = Tajima), span = 0.1)+
                        geom_vline(data = Signif_SNP %>%
                                        filter(scaf == signif_windows[signif_windows$Method == metric,]$chr[win] &
                                                pos >= signif_windows[signif_windows$Method == metric,]$Start[win] &
                                                pos <= signif_windows[signif_windows$Method == metric,]$End[win],
                                                color >0) %>%
                                        mutate(Y = max(sub_theta$Tajima)),
                                mapping = aes( xintercept=(pos/10^6), color = color), col = "red",size = 0.1, inherit.aes = FALSE)+
                        geom_hline(yintercept = 0, color = "black")+
                        theme_bw() +
                        labs(x = "Windows center (Mb)", y= "Tajima")


  list_plot[[paste(signif_windows[signif_windows$Method == metric,][win], collapse = "_")]] <- p1/p3 +
    plot_layout(heights = c(5,2))


}

Plot_taj <- wrap_plots(list_plot) + plot_layout(guide = "collect")
ggsave(paste("99_plots/Zoom_windows_lb_",metric,".png",sep =""),
       Plot_taj,
       device = "png",
       dpi = 300,
       width = 16.9, height = 11.38, units = "cm",)
ggsave(paste("99_plots/Zoom_windows_lb_",metric,".pdf",sep =""),
       Plot_taj,
       device = "pdf",
       dpi = 300,
       width = 16.9, height = 11.38, units = "cm",)
}
