#' plotCoverage
#'
#' Function to plot raw coverage for a tumour/normal sample
#' @param counts_file File to plot. [Required]
#' @param readLength Specifiy the length of read [Default 100]
#' @param windowSize Specifiy the window size used [Default 50000]
#' @import tidyr
#' @keywords plot chrom
#' @export
#' @examples plotCoverage(counts_file = 'data/counts/A785-A788R11.tagged.SC.hits.filt-vs-A785-A788R12.tagged.SC.hits.filt.window-50000.minw-4.count', readLength=100,windowSize = 50000)

plotCoverage <- function(counts_file = 'data/counts/A785-A788R11.tagged.SC.hits.filt-vs-A785-A788R12.tagged.SC.hits.filt.window-50000.minw-4.count',
                         rollmean=FALSE,
                         write=FALSE,
                         readLength=150,
                         chroms=c('2L', '2R', '3L', '3R', 'X'),
                         windowSize = 'auto'){
  counts <- read.delim(counts_file, header = T)

  if(windowSize=='auto') windowSize <- head(counts, 1)$end


  cat("Specified read length:", readLength , "\n")
  cat("Read counts are from a window size of:", windowSize , "\n")


  dir.create(file.path("plots"), showWarnings = FALSE)
  dir.create(file.path("plots", "coverage"), showWarnings = FALSE)

  base = basename(counts_file)
  parts <- strsplit(base, "[.]")[[1]]
  sample <- parts[1]
  plot_width <- 20
  plot_height <- 20

  covOut <- paste(sample, "coverage.png", sep='_')

  cat("Plotting coverage for ", sample, "\n", sep="")

  colours<-c( "#E7B800", "#00AFBB")

  df <- counts %>%
    gather(test, ref, key="sample", value = read_count) %>%
    dplyr::mutate(sample = recode(sample,
                           test = "Tumour",
                           ref = "Control")) %>%
    dplyr::mutate(depth = (read_count/windowSize)*readLength) %>%
    dplyr::mutate(pos = start/1e6) %>%
    dplyr::select(sample, chromosome, depth, pos) %>%

    dplyr::filter(chromosome %in% chroms) %>%
    dplyr::arrange(chromosome, pos)

  if(is.numeric(rollmean)){
    cat("Using a rolling mean every", rollmean , "windows\n")

    colours<-c( "#259FBF" )

    averageCounts <- df %>%
      filter(sample == 'Tumour') %>%
      mutate(moving_average = roll_mean(depth, rollmean, align="right", fill=0))

    p <- ggplot(averageCounts) +
      geom_density(aes(pos, moving_average, fill=sample, colour=sample), stat='identity', alpha =0.6) +
      scale_x_continuous("Mb", expand = c(0.001, 0.001), breaks = seq(0, max(averageCounts$pos), by = 1 )) +
      scale_y_continuous("Depth", expand = c(0.01, 0.01))
    p <- p + facet_wrap(~chromosome, ncol=1, scales = 'free_x')
    if(length(chroms) > 1) {
      p <- p + facet_wrap(~chromosome, ncol=1, scales = 'free_x')
    } else {
      plot_width <- 20
      plot_height <- 5
    }

    p <- p + cleanTheme() +
    theme( panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
           strip.text = element_text(size=20),
           axis.text = element_text(size=18),
           axis.title = element_text(size=20)) +
    scale_fill_manual(values=colours) +
    scale_colour_manual(values=colours)
  } else {
    p <- ggplot(df) +
      geom_density(aes(pos, depth, fill=sample, colour=sample), stat='identity', alpha =0.6) +
      scale_x_continuous("Mb", expand = c(0.001, 0.001), breaks = seq(0, max(df$pos), by = 1 )) +
      scale_y_continuous("Depth", expand = c(0.01, 0.01)) +

      facet_wrap(~chromosome, ncol=1, scales = 'free_x') +
      cleanTheme() +
        theme( panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
               strip.text = element_text(size=20),
               axis.text.y = element_text(size=15),
               axis.title = element_text(size=20)) +

      scale_fill_manual(values=colours) +
      scale_colour_manual(values=colours)
  }

  if(write){
    ggsave(paste("plots/coverage/", covOut, sep=""), width = plot_width, height = plot_height)
  } else {
    print(p)
  }
}
