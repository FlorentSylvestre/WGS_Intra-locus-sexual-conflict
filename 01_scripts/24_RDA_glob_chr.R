#Package
 library(data.table)
 library(ggplot2)
 library(magrittr)
 library(tidyr)
 library(dplyr)
 library(ggnewscale)
 chrconv <- data.frame("ARABE" = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21),
                       "ROMAIN" = c("chrI","chrII","chrIII","chrIV","chrV",
                                    "chrVI", "chrVII","chrVIII", "chrIX","chrX","chrXI","chrXII",
                                    "chrXIII","chrXIV","chrXV","chrXVI","chrXVII",
                                    "chrXVIII", "chrXX", "chrXXI"))
dat = fread("09_RDA/OUTPUT_RDA.pval")

dat %<>% filter(Scale != "PER_WIN")

for(i in 1:nrow(dat)){
   if(!is.na(dat$Pos[i])){ dat$Pos[i]<-  chrconv$ROMAIN[chrconv$ARABE == dat$Pos[i]]}
}

write.table(dat, "/99_plots/RDA_Result_table",
            quote = F,
            sep =  "\t",
            row.names = FALSE
            )
