#' allPlot
#'
#' Function to plot all chromosomes for all samples in path provided
#' @param path Sets the path to cnv files. [Default "data/"]
#' @import RColorBrewer
#' @import scales
#' @import ggplot2
#' @keywords plot all
#' @export
#' @examples allPlot(path="data/")


allPlot <- function(path = 'data/') {
  dir.create(file.path("plots"), showWarnings = FALSE)
  file.names <- dir(path, pattern = ".cnv")

  for (i in 1:length(file.names)) {
    cat("Processing file", file.names[i], "\n")
    parts <- strsplit(file.names[i], "[.]")[[1]]
    sample <- parts[1]

    read_file_in <- read.delim(paste("cnvs/", file.names[i], sep = ""), header = T)
    clean_file <- cleanR(read_file_in)
    cols <- brewer.pal(n = 7, name = "RdBu")

    ylim<-c(-5,5)

    p <- ggplot()
    p <- p + geom_point(data=clean_file, aes(start/1000000, log2, colour = log2), size = 1)
    p <- p + scale_alpha(range = c(0.1, 5))
    p <- p + scale_x_continuous("Mb", expand = c(0.001, 0.001))
    p <- p + scale_y_continuous("Log2 FC ratio",limits=ylim, expand = c(0, 0), breaks = seq(min(ylim),max(ylim),by=1))

    p <- p + scale_colour_gradientn(colours = cols, values = rescale(c(-3, -0.585, -0.001, 0, 0.001, 0.585, 3)), guide = "colorbar", limits = ylim)
    p <- p + facet_wrap(~chromosome, scale = "free_x", ncol = 2)
	
    p <- p + cleanTheme(base_size = 10)
    
    p <- p + ggtitle(paste(sample))
	
    outfile <- paste(sample, "_", "CNVs", ".pdf", sep = "")
    cat("Writing file", outfile, "to `plots/`", "\n")
    ggsave(paste("plots/", outfile, sep = ""), width = 20, height = 10)

    }
}