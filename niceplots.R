finaluno=cowplot::plot_grid(A,B,C,e = "AUTO",nrow=1,ncol=3,rel_widths = c(1,1,1),rel_heights=c(1,1,1),align="hv")

tiff("C:/Users/oskar/Downloads/figure2.tiff", units="px", width=1800*2, height=1476, res=300)
finaluno
dev.off()
