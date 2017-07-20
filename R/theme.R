#' clean.theme
#'
#' Clean theme for plotting
#' @param base_size = sets the base font size, defaults to 12
#' @import ggplot2
#' @keywords theme
#' @export
#' @examples
#' clean.theme()

clean_theme <- function(base_size = 12){
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text = element_text(size=12),
    axis.title = element_text(size=15)
    )
}