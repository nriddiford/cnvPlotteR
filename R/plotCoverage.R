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

plotCoverage <- function(counts_file = 'data/counts/A785-A788R11.tagged.SC.hits.filt-vs-A785-A788R12.tagged.SC.hits.filt.window-50000.minw-4.count', readLength=100,windowSize = 10000){
  counts <- read.delim(counts_file, header = T)

  cat("Specified read length:", readLength , "\n")

  dir.create(file.path("plots"), showWarnings = FALSE)
  dir.create(file.path("plots", "coverage"), showWarnings = FALSE)

  base = basename(counts_file)
  parts <- strsplit(base, "[.]")[[1]]
  sample <- parts[1]

  cat("Plotting coverage for ", sample, "\n", sep="")

  colours<-c( "#E7B800", "#00AFBB")

  df <- counts %>%
    gather(test, ref, key="sample", value = read_count) %>%
    dplyr::mutate(sample = recode(sample,
                           test = "Tumour",
                           ref = "Control")) %>%
    dplyr::mutate(depth = (read_count/windowSize)*readLength) %>%
    dplyr::mutate(pos = start/1000000) %>%
    dplyr::select(sample, chromosome, depth, pos) %>%

    dplyr::filter(chromosome != "Y", chromosome != 4) %>%
    dplyr::arrange(chromosome, pos)


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

    covOut<-paste(sample, "coverage.png", sep='_')
    cat("Writing file", covOut, "\n")
    ggsave(paste("plots/coverage/", covOut, sep=""), width = 20, height = 20)

    p

}
