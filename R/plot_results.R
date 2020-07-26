#------------------------------------------------
#' @title Plot sensitivity analysis
#'
#' @description The function is used to plot any sensitivity analysis that is
#' run. It takes in the output of \code{\link{coi_test}} and creates a grid of
#' plots.
#'
#' @param data The data to be plotted.
#' @param dims A list representing the number of rows and columns our plots
#' will be split into.
#' @param sub_title A list of titles for each individual subplot.
#' @param title The title of the overall figure.
#' @param caption The caption of the overall figure.
#'
#' @export

sensitivity_plot <- function(data, dims, sub_title = NULL, title = NULL,
                             caption = NULL){

  # Ensure that ggplot2 and ggpubr are installed
  if (!requireNamespace("ggplot2", quietly = TRUE) &
      !requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Packages \"ggplot2\" and \"ggpubr\" must be installed in order to plot the tests.",
         call. = FALSE)
  } else if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" must be installed in order to plot the tests.",
         call. = FALSE)
  } else if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package \"ggpubr\" must be installed in order to plot the tests.",
         call. = FALSE)
  }

  # Check inputs
  assert_eq(names(data),
            c("predicted_coi", "probability", "param_grid", "boot_error"))
  assert_pos_int(dims, zero_allowed = FALSE)
  assert_length(dims, 2)
  if (!is.null(sub_title)) {assert_vector(sub_title)}
  if (!is.null(title)) {assert_single_string(title)}
  if (!is.null(caption)) {assert_single_string(caption)}

  # Convert the predicted_coi dataframe into a long format. More specifically,
  # establish a column for the true coi and a column for the estimated coi
  # value. In some cases, we change parameters other than the COI. When this is
  # the case, we want to make multiple plots and see the effect of the changing
  # parameter. Thus, we establish another column: loop number, that tells us
  # how often the other parameter is changing. In the end, there will be three
  # columns: true_COI, estimated_COI, and loop_number.
  plot_df <- data$predicted_coi %>%
    tidyr::gather("true_COI", "estimated_COI") %>%
    tidyr::extract(.data$true_COI, c("true_COI", "loop_number"),
                   "coi_(.+)_(.+)") %>%
    dplyr::mutate_all(as.numeric)

  # We determine how many different panels there will be by finding the unique
  # values of loop_number
  num_loops <- unique(plot_df$loop_number)

  # Ensure that there are enough panels to include all the graphs
  user_dims = dims[1] * dims[2]
  needed_dims    = length(num_loops)
  if (!all(user_dims >= needed_dims)) {
    message <- sprintf("Not enough panels have been specified. User input %s
                       panel(s), but %s panel(s) needed.", user_dims,
                       needed_dims) %>%
      stringr::str_squish() %>%
      stringr::str_wrap()
    stop(message, call. = FALSE)
  }

  # We then call a helper function: sensitivity_plot_element, that creates each
  # individual plot and store these plots as a list
  myplots <- lapply(num_loops,
                    sensitivity_plot_element,
                    data = plot_df,
                    sub_title = sub_title)

  # Arrange the plots
  if (dims[1] * dims[2] == 1){
    # Do not include a panel label if there is only one plot
    arranged_plots <- ggpubr::ggarrange(plotlist = myplots,
                                        nrow = dims[1],
                                        ncol = dims[2])
  } else {
    # Include panel labels if there are more than one plots
    arranged_plots <- ggpubr::ggarrange(plotlist = myplots,
                                        labels = "AUTO",
                                        font.label = list(size = 10),
                                        nrow = dims[1],
                                        ncol = dims[2])
  }


  arranged_plots <-
    ggpubr::annotate_figure(arranged_plots,
                            top = ggpubr::text_grob(title, size = 13),
                            bottom = ggpubr::text_grob(caption,
                                                       hjust = 0,
                                                       x = 0.01,
                                                       size = 10))

  # Lastly, we return the arranged plots
  return(arranged_plots)
}

#------------------------------------------------
#' @title Plot a single panel of the sensitivity analysis
#'
#' @description Plots a single panel of the sensitivity analysis. Used as a
#' helper function to \code{\link{sensitivity_plot}}.
#'
#' @param loop_num The loop number. Represents how many total panels will
#' be plotted.
#' @inheritParams sensitivity_plot
#'
#' @keywords internal

sensitivity_plot_element <- function(data,
                                     loop_num,
                                     sub_title){

  # Ensure that ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop('Package \"ggplot2\" must be installed in order to plot the test.',
       call. = FALSE)
  }

  # Plot the figure
  single_plot <-
    ggplot2::ggplot(dplyr::filter(data, .data$loop_number == loop_num),
                    ggplot2::aes(x = .data$true_COI, y = .data$estimated_COI)) +
    ggplot2::geom_count(color = "blue", alpha = 0.7, show.legend = FALSE) +
    ggplot2::scale_size_area() +
    ggplot2::geom_abline(color = "red", size = 1) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title   = ggplot2::element_text(hjust = 0.5, size = 10),
                   axis.title   = ggplot2::element_text(size = 7),
                   legend.title = ggplot2::element_text(size = 7),
                   legend.text  = ggplot2::element_text(size = 7)) +
    ggplot2::labs(x = "True COI", y = "Estimated COI", title = sub_title[loop_num])

  return(single_plot)
}

#------------------------------------------------
#' @title Plot the mean absolute errors for a test
#'
#' @description Plots the mean absolute error and the confidence interval.
#'
#' @param data The data to be plotted.
#' @param fill The variable the data will be separated by.
#' @param fill_levels The levels for the fill variable.
#' @param title The title of the plot. Default to \code{NULL}.
#' @param legend_title The text for the legend. Default to \code{NULL}.
#' @param legend.position The position of the legend. One of \code{"none",
#' "left", "right", "bottom", "top"}.
#' @param second_fill Indicates if there will be a second fill variable and
#' what it will be.
#'
#' @export

error_plot <- function(data, fill = "COI", fill_levels = NULL, title = NULL,
                       legend_title = fill, legend.position = "bottom",
                       second_fill = NULL){

  # Ensure that ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" must be installed in order to plot the test.",
         call. = FALSE)
  }

  # Check inputs
  assert_eq(names(data),
            c("predicted_coi", "probability", "param_grid", "boot_error"))
  assert_single_string(fill)
  if (!is.null(fill_levels)) {assert_string(fill_levels)}
  if (!is.null(title)) {assert_single_string(title)}
  assert_single_string(legend_title)
  assert_single_string(legend.position)
  if (!is.null(second_fill)) {assert_single_string(second_fill)}

  # Convert data to a tibble
  plot_data <- data$boot_error %>% tidyr::unnest(cols = names(data$boot_error))

  # Make several changes to the data frame. For columns that have more than 1
  # unique value and are not one of c("mae", "lower", "upper", "bias"), convert
  # them to a factor. We also replace 0 with NA for the columns mae, lower, and
  # upper. This makes plotting easier and ensures that error bars are not shown
  # for data points where the mean absolute error is 0.
  plot_data <- plot_data %>%
    dplyr::mutate(dplyr::across(where(function(x) dplyr::n_distinct(x) > 1) &
                                  !tidyselect::all_of(c("mae", "lower",
                                                        "upper", "bias")),
                                as.factor)) %>%
    dplyr::mutate(mae   = dplyr::na_if(.data$mae, 0)) %>%
    dplyr::mutate(lower = dplyr::na_if(.data$lower, 0)) %>%
    dplyr::mutate(upper = dplyr::na_if(.data$upper, 0))

  # Customize the levels of the fill variable
  if (!is.null(fill_levels)){
    # Ensure that the number of levels input (fill_levels) are the same as the
    # number of levels for the fill variable
    if (length(fill_levels) != nlevels(plot_data[[fill]])){
      message <- sprintf("The length of the \"fill_levels\" variable (%s) does
                         not match with the number of levels of the fill
                         variable (%s)", length(fill_levels),
                         nlevels(plot_data[[fill]])) %>%
        stringr::str_squish() %>%
        stringr::str_wrap()
      stop(message, call. = FALSE)
    }
    else {
      # Customize labels
      levels(plot_data[[fill]]) <- fill_levels
    }
  }

  # Plot the data and return
  error_plot <- ggplot2::ggplot(plot_data,
                                ggplot2::aes(x = .data$COI, y = .data$mae,
                                             fill = eval(parse(text = fill)))) +
    ggplot2::geom_col(position = "dodge", na.rm = T) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                           width = .2, position = ggplot2::position_dodge(.9)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = legend.position,
                   plot.title   = ggplot2::element_text(hjust = 0.5, size = 13),
                   axis.title   = ggplot2::element_text(size = 10),
                   legend.title = ggplot2::element_text(size = 10),
                   legend.text  = ggplot2::element_text(size = 8)) +
    ggplot2::labs(y = "Mean Absolute Error", title = title, fill = legend_title)

  if (!is.null(second_fill)){
    error_plot <- error_plot + ggplot2::facet_wrap(~ plot_data[[second_fill]],
                                                   nrow = 2)
  }

  return(error_plot)

}

#------------------------------------------------
#' @title Plot a world map showing the COI
#'
#' @description Plot a world map showing the COI in each region where reads
#' were sampled from.
#'
#' @param data The data to be plotted.
#' @param variable The variable the data will plot.
#' @param label The label for the variable.
#' @param alpha The alpha value for the plotted data.
#' @param breaks The breaks for the color scale.
#'
#' @export
#'
world_map <- function(data, variable, label = NULL, alpha = 0.1,
                      breaks = c(1,2)){

  # Access world map data from ggplot2
  world <- ggplot2::map_data("world")

  # Plot world map
  map <- ggplot2::ggplot() +
    ggplot2::borders("world") +
    ggplot2::geom_polygon(data = world,
                          ggplot2::aes(x = .data$long, y = .data$lat,
                                       group = .data$group),
                          fill = "grey", alpha = 0.3) +
    ggplot2::geom_point(data = data,
                        ggplot2::aes(x = .data$long, y = .data$lat,
                                     size = variable, color = variable),
                        alpha = alpha) +
    ggplot2::scale_colour_viridis_c(limits = c(breaks[1],
                                               breaks[length(breaks)]),
                                    breaks = breaks) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_size(guide = "none") +
    ggplot2::labs(color = label) +
    ggplot2::coord_quickmap(xlim = c(-75, 150), ylim = c(-20, 20))

  return(map)
}
