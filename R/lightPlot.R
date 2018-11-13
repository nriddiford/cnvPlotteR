#' lightPlot
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
#' @import scales ggplot2 dplyr RColorBrewer
#' @keywords plot region
#' @export
#' @examples regionPlot(cnv_file="data/w500/test.window-500.cnv", from=3050000, to=3450000, chrom="X", ylim=c(-7,7), bp1=3129368,bp2=3352041, tick=100000, title="222Kb DEL on X")


lightPlot <- function(cnv_file, from=NULL, to=NULL, chrom=NULL, ylim=c(-5,5),
                      tick=100000, bp1=NULL, bp2=NULL, title=NULL, position=NULL) {

  if (missing(cnv_file) ) {
    if (!missing(position) || missing(from) || missing(to) || missing(chrom)){
      stop("Required arguments: `cnv_file` `from`, `to`, `chrom`", call. = FALSE)
    }
  }

  if (!missing(position)){
    vars <- unlist(strsplit(position, ":|-"))
    chrom <- vars[1]
    bp1 <- as.integer(vars[2])
    bp2 <- as.integer(vars[3])
    length <- bp2 - bp1
    from <- bp1 - length*5
    to <- bp2 + length*5
  }

  cat("Processing", cnv_file, "\n")
  cat("Specified region", from, "-", to, "on", chrom, "\n")

  if(to > 2500000 && from < 3750000 && chrom == "X"){
    cat("Specified Notch locus\n")
    notch <- TRUE
  } else {
    notch <- FALSE
  }

  if ( !missing(bp1) && !missing(bp2) || !missing(position) ){
    cat("Drawing lines for breakpoints: bp1:", bp1, " bp2:", bp2, "\n")
    draw_bps <- TRUE
  } else {
    draw_bps <- FALSE
  }

  base = basename(cnv_file)
  parts <- strsplit(base, "[.]")[[1]]
  sample <- parts[1]

  cat("Light-plotting region ", chrom, ":", from, "-", to, " for ", sample, "\n", sep="")

  read_file_in <- data.table::fread(cnv_file, header = T, select=c('chromosome', 'position', 'log2'))
  clean_file <- cleanR(read_file_in, region=T)

  if (missing(title)) title <- paste0(sample, chrom)

  # For white background
  customTheme = cleanTheme()
  barfill = 'black'

  from <- round_any(from, 1e5, floor)
  to <- round_any(to, 1e5, ceiling)

  region <- clean_file %>%
    dplyr::filter(chromosome == chrom) %>%
    dplyr::filter(position >= from,
                  position <= to) %>%
    mutate(mavlog2 = roll_mean(log2, 10, fill=0))

  ylim <- round_any(max(abs(region$mavlog2)), 1, ceiling)
  if (ylim < 2)
    ylim = 2
  ylim <- c(-ylim, ylim)

  p <- ggplot(region, aes(x = position, y = mavlog2))
  p <- p + geom_line(aes(colour = chromosome), size=2, show.legend = FALSE)
  p <- p + scale_x_continuous("Bp", expand = c(0.001, 0.001), breaks = seq(from, to, by = 50000))
  p <- p + scale_y_continuous("Log2 FC ratio", limits=c(min(ylim),max(ylim)), expand = c(0.02, 0.02), breaks = seq(min(ylim), max(ylim),by=1))
  p <- p + customTheme

  if(notch){
    kirreStart = 2740384
    dunceEnd = 3343000
    if (kirreStart < from){
      kirreStart = from
    }
    if (dunceEnd > to){
      dunceEnd = to
    }
    p <- p + annotate("rect", xmin=kirreStart, xmax=3134000, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.75, fill="#CFAEAEFE")
    p <- p + annotate("rect", xmin=3134000, xmax=3172000, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.75, fill="#8FBD80FE")
    p <- p + annotate("rect", xmin=3176000, xmax=dunceEnd, ymin=(min(ylim)+0.5), ymax=min(ylim), alpha=.75, fill="#A9D0DEFE")

    p <- p + annotate("text", x = kirreStart+10000, y = (min(ylim)+0.25), label="Kirre", size=6)
    p <- p + annotate("text", x = 3153000, y = (min(ylim)+0.25), label="Notch", size=6)
    p <- p + annotate("text", x = dunceEnd-100000, y = (min(ylim)+0.25), label="Dunce", size=6)
  }

  if(draw_bps){
    p <- p + geom_vline(xintercept = bp1, colour="slateblue", alpha=.7, linetype="dotted")
    p <- p + geom_vline(xintercept = bp2, colour="slateblue", alpha=.7, linetype="dotted")
  }

  scale_bar_start <- (from+(tick/10))
  scale_bar_end <- (scale_bar_start + tick)

  scale_text <- paste(tick/1000000, "Mb")
  p <- p + annotate("rect", xmin=scale_bar_start, xmax=scale_bar_end, ymin=(min(ylim)+0.6), ymax=(min(ylim)+0.7), fill=barfill)
  p <- p + annotate("text", x=(scale_bar_start + (tick/2)), y = (min(ylim)+0.9), label=scale_text, size=8, colour=barfill)

  p <- p + ggtitle(title)

  p
}
