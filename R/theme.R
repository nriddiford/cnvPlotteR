#' cleanTheme
#'
#' Clean theme for plotting
#' @param base_size Base font size [Default 12]
#' @import ggplot2
#' @keywords theme
#' @export

# Standard cleantheme for white backgrounds
cleanTheme <- function(base_size = 12){
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text = element_text(size=20),
    axis.title = element_text(size=20)
    )
}


# Modified for black background
blackTheme <- function(base_size = 12){
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, colour = 'white'),
    panel.background = element_blank(),
    plot.background = element_rect(colour = "black", fill = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="white", size = 0.5),
    axis.line.y = element_line(color="white", size = 0.5),
    axis.text = element_text(size=20, colour = 'white'),
    # axis.title = element_text(size=20, colour = 'white'),
    axis.ticks = element_line(color = "white"),
    axis.title = element_text(size = 20, color = "white"),
    legend.background = element_rect(color = NA, fill = "black"),
    legend.key = element_rect(color = "white",  fill = "black"),
    legend.text = element_text(size = base_size*0.8, color = "white"),
    legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white")
  )
}
