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
#' @import scales ggplot2 dplyr plyr RColorBrewer
#' @keywords plot region
#' @export
#' @examples regionPlot(cnv_file="data/w500/test.window-500.cnv", from=3050000, to=3450000, chrom="X", ylim=c(-7,7), bp1=3129368,bp2=3352041, tick=100000, title="222Kb DEL on X")


regionPlot <- function(cnv_file, from=2950000, to=3400000, chrom="X", ylim=c(-5,5),
                       tick=100000, bp1=NULL, bp2=NULL, position=NULL, title=NA, ext='png', theme = 'white', write=F) {

  if (missing(cnv_file) ) {
    if (!missing(position) || missing(from) || missing(to) || missing(chrom)){
      stop("Required arguments: `cnv_file` `from`, `to`, `chrom`", call. = FALSE)
    }
  }

  cat("Processing", cnv_file, "\n")

  if(!grepl('\\.', ext)) ext <- paste0('.', ext)
  cat("Specified region", from, "-", to, "on", chrom, "\n")

  if (!missing(position)){
    vars <- unlist(strsplit(position, ":|-"))
    chrom <- vars[1]
    bp1 <- as.integer(vars[2])
    bp2 <- as.integer(vars[3])
    length <- bp2 - bp1
    from <- bp1 - length*5
    to <- bp2 + length*5
  }

  if(from > 2500000 && from < 3172000 && to >= 3134000 && to < 3750000 && chrom == "X"){
    cat("Specified Notch locus\n")
    notch <- T
  } else {
    notch <- F
  }

  if ( !missing(bp1) && !missing(bp2) ){
    cat("Drawing lines for breakpoints: bp1=", bp1, " bp2=", bp2, "\n")
    draw_bps <- TRUE
  } else {
    draw_bps <- FALSE
  }

  cat("Chrom:", chrom, "\n")
  cat("Specified ylim", ylim , "\n")
  cat("Plotting ticks on", chrom, "every", tick, "\n")

  dir.create(file.path("plots"), showWarnings = FALSE)
  dir.create(file.path("plots", "regions"), showWarnings = FALSE)

  base = basename(cnv_file)
  parts <- strsplit(base, "[.]")[[1]]
  sample <- parts[1]

  cat("Plotting region ", chrom, ":", from, "-", to, " for ", sample, "\n", sep="")

  read_file_in <- read.delim(cnv_file, header = T)
  clean_file <- cleanR(read_file_in, region=T)

  if (missing(title)) title <- paste(sample, chrom, sep = " ")

  # For white background
  if(theme=='white'){
    cols <- c("#941212FE", "#C44747FE", "#B3A5A5FE", "#4FA9BDFE", "#248DB3FE")
    customTheme = cleanTheme()
    barfill = 'black'
  }
  # For black background
  if(theme=='black'){
    # cols <- c("#EB7609", "#EBA17C", "#B3A5A5FE", "#4FA9BDFE", "#248DB3FE")
    # cols <- c("#F2EC3F", "#EBDE96", "#B3A5A5FE", "#4FA9BDFE", "#248DB3FE")
    cols <- c("#F5D520", "#EBDC86", "#BCC1CC", "#9ACAD6", "#4CA8E6")
    customTheme = blackTheme()
    ext <- paste0('_dark', ext)
    barfill = 'white'
  }

  from <- plyr::round_any(from, 1000, floor)
  to <- plyr::round_any(to, 1000, ceiling)


  region <- clean_file %>%
    filter(chromosome == chrom) %>%
    filter(position >= from & position <= to) %>%
    droplevels()

  p <- ggplot(region, aes(start, log2, colour = log2))
  p <- p + geom_point(size=2.5)
  p <- p + scale_alpha(range = c(0.1, 5))
  p <- p + scale_x_continuous("Bp", expand = c(0.001, 0.001), breaks = seq(from, to, by = tick))
  p <- p + scale_y_continuous("Log2 FC ratio",limits=ylim, expand = c(0, 0), breaks = seq(min(ylim),max(ylim),by=1))
  p <- p + scale_colour_gradientn(colours = cols, values = rescale(c(-3, -0.25, 0, 0.25, 3)), guide = "colorbar", limits = ylim)
  p <- p + customTheme

  if(notch){
    kirreStart = 2740384
    if (kirreStart < from)
      kirreStart = from + 100
    dncEnd <- 3343000
    if (dncEnd > to)
      dncEnd = to - 100
    p <- p + annotate("rect", xmin=kirreStart, xmax=3134000, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.75, fill="#CFAEAEFE")
    p <- p + annotate("rect", xmin=3134000, xmax=3172000, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.75, fill="#8FBD80FE")
    p <- p + annotate("rect", xmin=3176000, xmax=dncEnd, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.75, fill="#A9D0DEFE")

    p <- p + annotate("text", x = kirreStart+(3134000-kirreStart)/2, y = (min(ylim)+0.25), label="Kirre", size=6)
    p <- p + annotate("text", x = 3153000, y = (min(ylim)+0.25), label="Notch", size=6)
    p <- p + annotate("text", x = 3176000+(dncEnd-3176000)/2, y = (min(ylim)+0.25), label="Dunce", size=6)
  }

  if(draw_bps){
    p <- p + geom_vline(xintercept = bp1, colour="slateblue", alpha=.7, linetype="dotted")
    p <- p + geom_vline(xintercept = bp2, colour="slateblue", alpha=.7, linetype="dotted")
  }

  scale_bar_start <- (from+(tick/10))
  scale_bar_end <- (scale_bar_start + tick/10)
  text_pos <- scale_bar_start + (scale_bar_end - scale_bar_start)


  scale_text <- paste(tick/1e6, "Mb")
  p <- p + annotate("rect", xmin=scale_bar_start, xmax=scale_bar_end, ymin=(min(ylim)+1), ymax=(min(ylim)+0.8), fill=barfill)
  p <- p + annotate("text", x=text_pos, y = (min(ylim)+1.2), label=scale_text, size=8, colour=barfill)

  p <- p + ggtitle(title)

  if(notch){
    outfile <- paste0(sample, ".Notch", ext)
  } else {
    outfile <- paste0(sample, ".", chrom, "_", from, "-", to, ".", ext)
  }
  if(write){
    cat("Writing file", outfile, "to '../plots/regions/'", "\n")
    ggsave(paste("plots/regions/", outfile, sep = ""), width = 20, height = 10)
  }
  p
}
