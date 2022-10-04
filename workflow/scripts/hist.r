library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

# set args
args = commandArgs(trailingOnly=TRUE)
input <- args[1]
output_csv <- args[2]
output_png <- args[3]

#input <- "~/../../../Volumes/ccbr1214/v1.0/results/03_aligned/02_bed/Sample1.mapped.bed"
#output_csv <- "~/../../../Volumes/ccbr1214/v1.0/results/03_aligned/03_histograms/Sample1.length_hist_all.csv"
#output_png <- "~/../../../Volumes/ccbr1214/v1.0/results/03_aligned/03_histograms/Sample1.length_hist_all.png"

# read in rawdatam subset
rawdata <- read.csv(input, header=FALSE, sep = '\t')
data <- rawdata[,1:3]
colnames(data) <- c('Chrom', 'Start', 'End')

#calculate length
data$Length <- data$End - data$Start

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
ggsave(output_png,p1)

# pull data, save raw data
pg <- ggplot_build(p1)
write.csv(pg$data, output_csv)