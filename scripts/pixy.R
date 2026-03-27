library(ggplot2)
library(tidyverse)

# dxy plot function
plot_dxy <- function(file_name, my_colors, output_name) { 
  inp <- read.table(file_name, sep="\t", header=T) # read the data
  # we reorder the chromosome factor levels to ensure they are plotted in the correct order
  chroms <- unique(inp$chromosome)
  chrOrder <- sort(chroms)
  inp$chrOrder <- factor(inp$chromosome, levels=chrOrder) 
  inp$comparison <- paste(inp$pop1, "vs", inp$pop2) # create a new column for the population comparison
  # if the avg_dxy column is present, we proceed to plot
  if("avg_dxy" %in% colnames(inp)){
    pops <- unique(inp[c("pop1", "pop2")]) # we get the unique population names to assign colors

       # we plot dxy across the genome, faceting by population comparison and chromosome, and coloring by comparison
        popPlot <- ggplot(inp, aes(window_pos_1, y=avg_dxy, color=comparison)) +
            geom_point()+
            geom_smooth(color="black", method ="loess", se = FALSE, span=0.1)+ # we add a smoothed line to show the overall trend
            scale_x_continuous(labels = function(x) x / 1000000) +
            facet_grid(comparison ~ chrOrder, scales = "free_x", space = "free_x")+ 
            labs(title = paste("Dxy Comparison:", output_name), x = "Position (Mb)", y = "Dxy") +
            theme_classic(base_size = 18) + 
            theme(legend.position = "none", strip.text = element_text(face = "bold", size = 20),
              axis.title = element_text(face = "bold"),axis.text = element_text(size = 14),plot.title = element_text(face = "bold", size = 22, hjust = 0.5) )+
            scale_color_manual(values = my_colors)
            theme_classic()+
            theme(legend.position = "none")
        ggsave(paste("dxyplot_",output_name,".png", sep=""), plot = popPlot, device = "png",width = 16,height = 10, dpi = 300)
    }
else {
  print("Dxy not found in this file")}
}

# we call the function for each of our dxy files, specifying the colors and output names for the plots

# GQ>=10, window size 10kb
plot_dxy("results_10_10_dxy.txt",c("deeppink4","hotpink1","lightpink1"),"window_size10kb_gq10")

# GQ>=10, window size 50kb
plot_dxy("results_50_10_dxy.txt",c("orchid4","mediumorchid1","plum"),"window_size50kb_gq10")

# GQ>=15, window size 50kb
plot_dxy("results_50_15_dxy.txt",c("darkgreen","seagreen4","palegreen3"),"window_size50kb_gq15")

# GQ>=15, window size 10kb
plot_dxy("results_10_15_dxy.txt",c("royalblue4","lightslateblue","steelblue1"),"window_size10kb_gq15")

# GQ>=15, window size 100kb
plot_dxy("results_100_15_dxy.txt",c("royalblue4","lightslateblue","steelblue1"),"window_size100kb_gq15")
