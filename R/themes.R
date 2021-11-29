#------------------------------------------------
#' Custom ggplot2 theme
#'
#' @inheritParams ggplot2::theme_classic
#'
#' @importFrom ggplot2 %+replace%
#' @export
#' @examples
#' library("ggplot2")
#' p <- ggplot(mtcars, aes(x = wt, y = mpg, colour = factor(gear))) +
#'   geom_point() +
#'   facet_wrap(~am) +
#'   geom_smooth(method = "lm", se = FALSE)
#'
#' p + theme_coiaf()
theme_coiaf <- function(base_size = 10,
                        base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
      legend.title = ggplot2::element_text(size = 12),
      legend.position = "right",
      complete = TRUE
    )
}

# Need to have a function that returns the default theme which will be used in
# figures. This is according to the best practices for using ggplot2 in a
# package. Link:
# <https://ggplot2.tidyverse.org/articles/ggplot2-in-packages.html>
default_theme <- function() {
  theme_coiaf()
}
