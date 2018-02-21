#' regionPlot
#'
#' Function to plot a particular region from a single chromosomes for a sample
#' @param cnv_file File to plot. [Required]
#' @param from Start coordinate of plotting window [Default from=2950000 ]
#' @param to End coordinate of plotting window [Default to=3400000]
#' @param chrom Specify the chromosome to plot [Default 'X']
#' @param ylim Adjust y axis limits in plot [Default ylim=c(-5,5)]
#' @param tick Adjust spacing of ticks [Defualt 1Mb]
#' @param bp1 Draw a vertical line at given coordinate to mark bp1
#' @param bp2 Draw a vertical line at given coordinate to mark bp2
#' @param title Plot title
#' @import RColorBrewer
#' @import scales
#' @import ggplot2
#' @keywords plot region
#' @export
#' @examples regionPlot(cnv_file="data/w500/test.window-500.cnv", from=3050000, to=3450000, chrom="X", ylim=c(-7,7), bp1=3129368,bp2=3352041, tick=100000, title="222Kb DEL on X")


regionPlot <- function(cnv_file, from=NA, to=NA, chrom="X", ylim=c(-5,5), tick=100000, bp1=NA, bp2=NA, title=NA) {

  cat("Processing", cnv_file, "\n")

  if (is.na(from & to)) {
    from<-2950000
    to<-3400000
    cat("No region specified - defaulting to X:2950000-3400000", "\n")
    notch <- 1
  }
  else {
    cat("Specified region", from, "-", to, "on", chrom, "\n")
    notch <- F
  }

  if (!is.na(tick)) {
    cat("Plotting ticks on", chrom, "every", tick, "\n")
  }

  if (!is.na(bp1 & bp2)){
    cat("Drawing lines for breakpoints: bp1=", bp1, " bp2=", bp2, "\n")
    draw_bps <- 1
  }
  else {
    draw_bps <- F
  }


  cat("Chrom:", chrom, "\n")
  cat("Specified ylim", ylim , "\n")

  dir.create(file.path("plots"), showWarnings = FALSE)
  dir.create(file.path("plots", "regions"), showWarnings = FALSE)

  base = basename(cnv_file)
  parts <- strsplit(base, "[.]")[[1]]
  sample <- parts[1]

  cat("Plotting region ", chrom, ":", from, "-", to, " for ", sample, "\n", sep="")

  read_file_in <- read.delim(cnv_file, header = T)
  clean_file <- cleanR(read_file_in, region=T)
  # Change cols to give darker cols - 2.2.18
  cols <- c("#941212FE", "#C44747FE", "#B3A5A5FE", "#4FA9BDFE", "#248DB3FE")
  # cols <- brewer.pal(n = 5, name = "RdBu")

  chrom_data <- subset(clean_file, clean_file$chromosome == chrom)

  region <- subset(chrom_data, position>=from & position<=to)
  p <- ggplot(data=region, aes(start, log2, colour = log2))
  p <- p + geom_point(size=2.5)
  p <- p + scale_alpha(range = c(0.1, 5))
  p <- p + scale_x_continuous("Bp", expand = c(0.001, 0.001), breaks = seq(from, to, by = tick))
  p <- p + scale_y_continuous("Log2 FC ratio",limits=ylim, expand = c(0, 0), breaks = seq(min(ylim),max(ylim),by=1))
  #p <- p + scale_colour_gradientn(colours = cols, values = rescale(c(-3, -0.585, -0.001, 0, 0.001, 0.585, 3)), guide = "colorbar", limits = ylim)

  p <- p + scale_colour_gradientn(colours = cols, values = rescale(c(-3, -0.25, 0, 0.25, 3)), guide = "colorbar", limits = ylim)

  # p <- p + blackTheme()
  p <- p + cleanTheme()

  # Draw lines for log values corresponding to FC 2, 1.5 and 1.25
  p <- p + geom_hline(yintercept = -1, colour = "royalblue", alpha = 0.4)
  p <- p + geom_hline(yintercept = 1, colour = "royalblue", alpha = 0.4)
  p <- p + geom_hline(yintercept = -0.585, colour = "blue", alpha = 0.4)
  p <- p + geom_hline(yintercept = 0.585, colour = "blue", alpha = 0.4)
  p <- p + geom_hline(yintercept = -0.322, colour = "slateblue", alpha = 0.4, linetype = "dotted")
  p <- p + geom_hline(yintercept = 0.322, colour = "slateblue", alpha = 0.4, linetype = "dotted")

  if(notch){
    # p <- p + annotate("rect", xmin=2950000, xmax=3134000, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.075, fill="skyblue4")
    # p <- p + annotate("rect", xmin=3134000, xmax=3172000, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.075, fill="skyblue")
    # p <- p + annotate("rect", xmin=3176000, xmax=3343000, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.075, fill="slateblue")
    #kirre = 2740384
    p <- p + annotate("rect", xmin=2950000, xmax=3134000, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.75, fill="#CFAEAEFE")
    p <- p + annotate("rect", xmin=3134000, xmax=3172000, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.75, fill="#8FBD80FE")
    p <- p + annotate("rect", xmin=3176000, xmax=3343000, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.75, fill="#A9D0DEFE")

    p <- p + annotate("text", x = 3037000, y = (min(ylim)+0.25), label="Kirre", size=6)
    p <- p + annotate("text", x = 3153000, y = (min(ylim)+0.25), label="Notch", size=6)
    p <- p + annotate("text", x = 3259000, y = (min(ylim)+0.25), label="Dunce", size=6)
  }

  if(draw_bps){
    p <- p + geom_vline(xintercept = bp1, colour="slateblue", alpha=.7, linetype="dotted")
    p <- p + geom_vline(xintercept = bp2, colour="slateblue", alpha=.7, linetype="dotted")
  }


  scale_bar_start <- (from+(tick/10))
  scale_bar_end <- (scale_bar_start + tick)

  scale_text <- paste(tick/1000000, "Mb")
  p <- p + annotate("rect", xmin=scale_bar_start, xmax=scale_bar_end, ymin=(min(ylim)+1), ymax=(min(ylim)+0.8), fill="black")
  p <- p + annotate("text", x=(scale_bar_start + (tick/2)), y = (min(ylim)+1.2), label=scale_text, size=8)
  if (!is.na(title)){
    p <- p + ggtitle(paste(title))

  }
  else {
    p <- p + ggtitle(paste(sample, " ", chrom, sep = ""))
  }

  if(notch){
    outfile <- paste(sample, ".", "Notch", ".png", sep = "")
  }
  else {
    outfile <- paste(sample, ".", chrom, "_", from, "-", to, ".png", sep = "")
  }
  cat("Writing file", outfile, "to '../plots/regions/'", "\n")
  ggsave(paste("plots/regions/", outfile, sep = ""), width = 20, height = 10)
  p
}
