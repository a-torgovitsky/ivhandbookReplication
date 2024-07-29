# ggplot2 theme
base_theme <- function() {
  theme_bw() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal")
}
