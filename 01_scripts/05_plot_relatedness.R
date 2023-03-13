library(ggplot2)
library(magrittr)

data <- fread("05_relatedness/NO_INV_COV_MISS_QUAL_MAF010_filtered_99_autosomes_biallelic_SOR_sorted_renamed_NoUn.relatedness2")

heat <- data %>%
  ggplot(aes(x = INDV1, y =INDV2)) +
  geom_tile(aes(fill = RELATEDNESS_PHI)) +
  scale_fill_gradientn(colours=c("white", "lightgrey", "turquoise", "deepskyblue3", "blue", "blue3", "navyblue", "black"),
                       values=c(0, 0.125, 0.25, 0.375, 0.5), name = "Relatedness") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(x = "Individuals", y = "Individuals")

rect_h <-ggplot()+
    geom_rect(aes(xmin = 0, xmax = 50.5,ymin = 0, ymax = 1), fill = "red")+
  geom_rect(aes(xmin = 50.5, xmax = 99,ymin = 0, ymax = 1), fill = "blue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

rect_v <-ggplot()+
  geom_rect(aes(ymin = 0, ymax = 50.5, xmin = 0, xmax = 1), fill = "red")+
  geom_rect(aes(ymin = 50.5, ymax = 99,xmin = 0, xmax = 1), fill = "blue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())


library(patchwork)
layout = "
ABBB
#CCC
"
heatmap <- (rect_v + heat + rect_h  + plot_layout(design = layout, widths = c(0.1, 6),heights = c(10, 0.1)))

ggsave(file = "99_plots/Heatmap_KING_inference.png",
       device = "png",
       plot = heatmap,width = unit(10, "cm"),height = unit(10, "cm"),
       dpi = 300)

ggsave(file = "99_plots/Heatmap_KING_inference.pdf",
       device = "pdf",
       plot = heatmap,width = unit(10,"cm"),height = unit(10,"cm"),
       dpi = 300 )
