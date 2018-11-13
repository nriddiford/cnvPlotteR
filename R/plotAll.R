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


allPlot <- function(path = 'data/', outdir = 'plots') {
  dir.create(file.path(outdir), showWarnings = FALSE)
  file.names <- dir(path, pattern = ".cnv")

  for (i in 1:length(file.names)) {
    cat("Processing file", file.names[i], "\n")
    parts <- strsplit(file.names[i], "[.]")[[1]]
    sample <- parts[1]

    read_file_in <- read.delim(paste(path, file.names[i], sep = ""), header = T)
    clean_file <- cleanR(read_file_in)
    clean_file <- filter(clean_file, chromosome != "Y" & chromosome != "4")

    # cols <- brewer.pal(n = 7, name = "RdBu")
    cols <- c("#941212FE", "#C44747FE", "#B3A5A5FE", "#4FA9BDFE", "#248DB3FE")
    ylim<-c(-5,5)

    p <- ggplot()
    p <- p + geom_point(data=clean_file, aes(start/1000000, log2, colour = log2), size = 1)
    p <- p + scale_alpha(range = c(0.1, 5))
    p <- p + scale_x_continuous("Mb", expand = c(0.001, 0.001))
    p <- p + scale_y_continuous("Log2 FC ratio",limits=ylim, expand = c(0, 0), breaks = seq(min(ylim),max(ylim),by=1))

    # p <- p + scale_colour_gradientn(colours = cols, values = rescale(c(-3, -0.585, -0.001, 0, 0.001, 0.585, 3)), guide = "colorbar", limits = ylim)
    p <- p + scale_colour_gradientn(colours = cols, values = rescale(c(-3, -0.25, 0, 0.25, 3)), guide = "colorbar", limits = ylim)

    p <- p + facet_wrap(~chromosome, scale = "free_x", ncol = 2)

    # p <- p + blackTheme() +
    p <- p + cleanTheme() +
      theme(
        axis.text = element_text(size=15)
      )

    p <- p + ggtitle(paste(sample))

    outfile <- paste(sample, "_", "CNVs", ".png", sep = "")
    cat("Writing file", outfile, "to", outdir, "\n")
    ggsave(paste(outdir, outfile, sep = "/"), width = 20, height = 10)

    }
}
