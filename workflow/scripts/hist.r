library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]

rawdata <- read.csv(input, header=FALSE, sep = '\t')
data <- rawdata[,1:3]

colnames(data) <- c('Chrom', 'Start', 'End')
data$Length <- data$End - data$Start

p1 <- ggplot() +


        geom_histogram(data = data, aes(x=Length), color = 'black', fill = 'white', binwidth = 1, size = 1) +

        theme_classic(base_size = 20) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
            text = element_text(colour = "black"),
            rect = element_rect(fill = "transparent"),
            panel.grid.major = element_line(colour = "grey", size = 0.5),
            panel.grid.minor = element_line(colour = "transparent", size = 0.5),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA)) +

        labs(y='Counts', x='Length')


pg <- ggplot_build(p1)
write.csv(pg$data, output)