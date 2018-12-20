#' Function to plot a single chromosomes for a given sample
#' @param chrom Specify the chromosome to plot [Defualt 'X']
#' @param cnv_file File to plot. [Required]
#' @param ylim Adjust y axis limits in plot [Default ylim=c(-5,5)]
#' @import RColorBrewer
#' @import scales
#' @import ggplot2
#' @keywords plot chrom
#' @export
#' @examples chromPlot(cnv_file = "data/test.window-10000.cnv", chrom = "3R", ylim=c(-6,6))
chromPlot <- function(cnv_file, chrom = NA, ylim=c(-5,5), ext='png',
                      theme = 'white', write=TRUE) {

  cat("Processing", cnv_file, "\n")

  if (is.na(chrom)) {
    chrom <- "X"
    cat("No chromosome specified - defaulting to", chrom, "\n")
  }
  else {
    cat("Plotting", chrom, "chromosome", "\n")
  }

  dir.create(file.path("plots"), showWarnings = FALSE)
  dir.create(file.path("plots", "chroms"), showWarnings = FALSE)

  base = basename(cnv_file)
  parts <- strsplit(base, "[.]")[[1]]
  sample <- parts[1]

  read_file_in <- read.delim(cnv_file, header = T)
  clean_file <- cleanR(read_file_in)

  # For white background
  if(theme=='white'){
    cols <- c("#941212FE", "#C44747FE", "#B3A5A5FE", "#4FA9BDFE", "#248DB3FE")
    customTheme = cleanTheme()
  } else if(theme=='black'){ # For black background
    # cols <- c("#EB7609", "#EBA17C", "#B3A5A5FE", "#4FA9BDFE", "#248DB3FE")
    cols <- c("#F5D520", "#EBDC86", "#BCC1CC", "#9ACAD6", "#4CA8E6")
    customTheme = blackTheme()
    ext <- paste0('_dark.', ext)
  } else{
    cat(theme, " not a valid theme. Please use 'white' or 'black'", "\n")
    break
  }
  chrom_data <- dplyr::filter(clean_file, chromosome == chrom)

  p <- ggplot(chrom_data, aes(start/1000000, log2, colour = log2))
  p <- p + geom_point()
  p <- p + scale_alpha(range = c(0.1, 5))
  p <- p + scale_x_continuous("Mb", expand = c(0.001, 0.001), breaks = seq(0, (max(chrom_data$end)-1), by = 1))
  p <- p + scale_y_continuous("Log2 FC ratio",limits=ylim, expand = c(0, 0), breaks = seq(min(ylim),max(ylim),by=1))
  # p <- p + scale_colour_gradientn(colours = cols, values = rescale(c(-3, -0.585, -0.001, 0, 0.001, 0.585, 3)), guide = "colorbar", limits = ylim)
  p <- p + scale_colour_gradientn(colours = cols, values = rescale(c(-3, -0.25, 0, 0.25, 3)), guide = "colorbar", limits = ylim)

  p <- p + customTheme

  if(chrom == "X"){
    p <- p + geom_vline(xintercept = 3.134, colour="slateblue", alpha=.7, linetype="dotted") + geom_vline(xintercept = 3.172, colour="slateblue", alpha=.7, linetype="dotted")
  }

  p <- p + ggtitle(paste(sample, " ", chrom, sep = ""))

  if(write){
    outfile <- paste(sample, "_", chrom, "_", "CNVs.", ext, sep = "")
    cat("Writing file", outfile, "to '../plots/chroms/'", "\n")
    ggsave(paste("plots/chroms/", outfile, sep = ""), width = 20, height = 10)
  }
  p
}
