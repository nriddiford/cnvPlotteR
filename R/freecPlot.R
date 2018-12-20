#' FreecPlot
#'
#' Function to plot the corrected CN for Control Freec "ratio" files.
#' @param cnv_file File to plot. [Required]
#' @param chroms Chromosomes to plot
#' @param write Write plot to file? T/F
#' @import RColorBrewer
#' @import scales
#' @import dplyr
#' @keywords plot normalised CN for Control Freec ratio file
#' @export
#' @examples freecPlot(cnv_file = "data/freec/HUM-7.tagged.filt.SC.RG.bam_ratio.txt")

freecPlot <- function(cnv_file = 'data/freec/HUM-7.tagged.filt.SC.RG.bam_ratio.txt',
                      chroms=c('2L', '2R', '3L', '3R', 'X'),
                      write=F) {

  cat("Processing", cnv_file, "\n")

  base = basename(cnv_file)
  parts <- strsplit(base, "[.]")[[1]]
  sample <- parts[1]
  # chrom_lengths = read.delim(chrom_lengths_file, header = F)
  # colnames(chrom_lengths) <- c("chr", "length")

  ratio <- read.delim(cnv_file, header = T)
  sex_chroms <- grep(pattern = "X|Y", chroms, value = T)
  autosomes <- grep(pattern = "X|Y", chroms, value = T, invert = T)

  green <- '#12A116FE'
  sub_green <- '#9DD486FE'
  grey <- '#7A8B8B'
  red <- '#CC2121FE'
  sub_red <- '#DBD828FE'

  ratio <- ratio %>%
    dplyr::filter(Chromosome %in% chroms) %>%
    droplevels()

  ratio$colour <- ifelse(ratio$Chromosome %in% autosomes & ratio$Subclone_CN>2, sub_green, grey)
  ratio$colour <- ifelse(ratio$Chromosome %in% sex_chroms & ratio$Subclone_CN>1, sub_green, ratio$colour)

  ratio$colour <- ifelse(ratio$Chromosome %in% autosomes & ratio$Subclone_CN<2 & ratio$Subclone_CN > 0,  sub_red, ratio$colour)
  ratio$colour <- ifelse(ratio$Chromosome %in% sex_chroms & ratio$Subclone_CN<1 & ratio$Subclone_CN > 0, sub_red, ratio$colour)

  ratio$colour <- ifelse(ratio$Chromosome %in% autosomes & ratio$CopyNumber>2, green, ratio$colour)
  ratio$colour <- ifelse(ratio$Chromosome %in% sex_chroms & ratio$CopyNumber>1, green, ratio$colour)

  ratio$colour <- ifelse(ratio$Chromosome %in% autosomes & ratio$CopyNumber<2, red, ratio$colour)
  ratio$colour <- ifelse(ratio$Chromosome %in% sex_chroms & ratio$CopyNumber<1, red, ratio$colour)

  ratio$point_size <- ifelse(ratio$colour != grey, 0.9, 0.5)

  ploidy<-2
  p <- ggplot(ratio)
  p <- p + geom_point(aes(Start, Ratio*ploidy, colour=colour),stat='identity', alpha =0.5, size = ratio$point_size)
  p <- p + scale_y_continuous("Normalised copy number", limits = c(0,6), breaks = seq(0,6, by=1))
  p <- p + facet_wrap(~Chromosome, ncol=2, scales='free_x')
  p <- p + cleanTheme() +
    theme( panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
           strip.text = element_text(size=20),
           axis.text.y = element_text(size=15),
           axis.title = element_text(size=20))

  p <- p + scale_color_identity()

  if(write){
    dir.create(file.path("plots", "freec"), showWarnings = FALSE)
    freecOut <- paste(sample, "freec_cn.png", sep='_')
    cat("Writing file", freecOut, "\n")
    ggsave(paste("plots/freec/", freecOut, sep=""), width = 20, height = 10)
    }
  p

}
