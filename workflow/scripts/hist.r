library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

# set args
args = commandArgs(trailingOnly=TRUE)
input <- args[1]
output_csv <- args[2]

#input <- "~/../../../Volumes/ccbr1214/analysis/7_ARPE_5u_MNase/results/03_aligned/02_bed/7_ARPE_5u_MNase.hg19.120-140.2500_ALUdepleted.bed"
#output_csv <- "~/../../../Volumes/ccbr1214/analysis/7_ARPE_5u_MNase/results/03_aligned/03_histograms/7_ARPE_5u_MNase.hg19.length_hist.120-140.2500_ALUdepleted.csv"
#output_png <- "~/../../../Volumes/ccbr1214/analysis/7_ARPE_5u_MNase/results/03_aligned/03_histograms/7_ARPE_5u_MNase.hg19.length_hist.120-140.2500_ALUdepleted.png"
#input="~/../../../Volumes/Zhurkin-20/analysis/8_NB26_5_200/results/tmp/histo_frags/tmp2.bed"
#output_csv="~/../../../Volumes/Zhurkin-20/analysis/8_NB26_5_200/results/03_aligned/03_histograms/8_NB26_5_200.hg19.length_hist.all.csv"
#output_png="~/../../../Volumes/Zhurkin-20/analysis/8_NB26_5_200/results/03_aligned/03_histograms/8_NB26_5_200.hg19.length_hist.all.png"

# read in rawdata subset
data <- read.csv(input, header=FALSE, sep = ' ')
colnames(data) <- c('Chrom', 'Start', 'End', 'Length')

# create histogram, save
p1 <- ggplot() + geom_histogram(data = data, aes(x=Length), color = 'black', fill = 'white', binwidth = 1, size = 1) +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5), text = element_text(colour = "black", face="bold"), 
            rect = element_rect(fill = "transparent"),
            panel.grid.major = element_line(colour = "grey", size = 0.5),
            panel.grid.minor = element_line(colour = "transparent", size = 0.5),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA)) +
    labs(y='Counts', x='Length')
#ggsave(output_png,p1)

# pull data, save raw data
pg <- ggplot_build(p1)
write.csv(pg$data, output_csv)
