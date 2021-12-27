#------------------------------------------------
#' Sensitivity plot
#'
#' Creates a plot of the sensitivity analysis.
#'
#' \loadmathjax
#' The function takes in the output of [sensitivity()] and creates a grid of
#' plots. Each plot is created using [ggplot2::geom_count()]. The number of
#' observations at each location is counted and then the count is mapped to
#' point area on the plot.
#'
#' The x-axis is the true COI, and the y-axis is the estimated COI. The counts
#' are plotted in blue, and red line is drawn with the equation \mjseqn{y = x}.
#' This line indicates where the blue circles should be if the algorithm was
#' 100% correct.
#'
#' @param data The data to be plotted.
#' @param dims A list representing the number of rows and columns our plots
#' will be split into.
#' @param result_type An indicator that indicates if a count or boxplot should
#' be plotted.
#' @param sub_title A list of titles for each individual subplot.
#' @param title The title of the overall figure.
#' @param caption The caption of the overall figure.
#'
#' @seealso [ggplot2::geom_count()] for more information on count plots and the
#' [ggplot2 website](https://ggplot2.tidyverse.org/index.html).
#' @family plotting functions
#' @export
sensitivity_plot <- function(data,
                             dims,
                             result_type,
                             sub_title = NULL,
                             title = NULL,
                             caption = NULL) {

  # Ensure that ggplot2 and ggpubr are installed
  if (!requireNamespace("ggplot2", quietly = TRUE) &
    !requireNamespace("ggpubr", quietly = TRUE)) {
    stop(
      "Packages \"ggplot2\" and \"ggpubr\" must be installed in order to plot the tests.",
      call. = FALSE
    )
  } else if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed in order to plot the tests.",
      call. = FALSE
    )
  } else if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop(
      "Package \"ggpubr\" must be installed in order to plot the tests.",
      call. = FALSE
    )
  }

  # Check inputs
  assert_in(
    names(data),
    c("predicted_coi", "probability", "param_grid", "boot_error")
  )
  assert_single_string(result_type)
  assert_in(result_type, c("disc", "cont"))
  assert_pos_int(dims, zero_allowed = FALSE)
  assert_length(dims, 2)
  if (!is.null(sub_title)) assert_vector(sub_title)
  if (!is.null(title)) assert_single_string(title)
  if (!is.null(caption)) assert_single_string(caption)

  # Convert the predicted_coi dataframe into a long format. More specifically,
  # establish a column for the true coi and a column for the estimated coi
  # value. In some cases, we change parameters other than the COI. When this is
  # the case, we want to make multiple plots and see the effect of the changing
  # parameter. Thus, we establish another column: loop number, that tells us
  # how often the other parameter is changing. In the end, there will be three
  # columns: true_coi, estimated_coi, and loop_number.
  plot_df <- data$predicted_coi %>%
    tidyr::gather("true_coi", "estimated_coi") %>%
    tidyr::extract(
      .data$true_coi, c("true_coi", "loop_number"),
      "coi_(.+)_(.+)"
    ) %>%
    dplyr::mutate(dplyr::across(.fns = as.numeric))

  # We determine how many different panels there will be by finding the unique
  # values of loop_number
  num_loops <- unique(plot_df$loop_number)

  # Ensure that there are enough panels to include all the graphs
  user_dims <- dims[1] * dims[2]
  needed_dims <- length(num_loops)
  if (!all(user_dims >= needed_dims)) {
    message <- glue::glue(
      "Must specify enough plotting panels:",
      "\n\u2139 {needed_dims} panels are required.",
      "\n\u2716 User specified {user_dims} panels."
    )
    stop(message, call. = FALSE)
  }

  # We then call a helper function: sensitivity_plot_element, that creates each
  # individual plot and store these plots as a list
  myplots <- lapply(
    num_loops,
    sensitivity_plot_element,
    data = plot_df,
    result_type = result_type,
    sub_title = sub_title
  )

  # Arrange the plots
  if (dims[1] * dims[2] == 1) {
    # Do not include a panel label if there is only one plot
    arranged_plots <- ggpubr::ggarrange(
      plotlist = myplots,
      nrow = dims[1],
      ncol = dims[2]
    )
  } else {
    # Include panel labels if there are more than one plots
    arranged_plots <- ggpubr::ggarrange(
      plotlist = myplots,
      labels = "AUTO",
      font.label = list(size = 10),
      nrow = dims[1],
      ncol = dims[2]
    )
  }

  ggpubr::annotate_figure(
    arranged_plots,
    top = ggpubr::text_grob(title, size = 13),
    bottom = ggpubr::text_grob(caption, hjust = 0, x = 0.01, size = 10)
  )
}

#------------------------------------------------
#' Single sensitivity plot
#'
#' Creates a single plot of the sensitivity analysis. Used as a helper function
#' to [sensitivity_plot()].
#'
#' @param loop_num The loop number. Represents how many total panels will
#' be plotted.
#' @inheritParams sensitivity_plot
#'
#' @keywords internal
sensitivity_plot_element <- function(data, loop_num, result_type, sub_title) {

  # Ensure that ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      'Package \"ggplot2\" must be installed in order to plot the test.',
      call. = FALSE
    )
  }

  # Plot the figure
  single_plot <- ggplot2::ggplot(
    dplyr::filter(data, .data$loop_number == loop_num),
    ggplot2::aes(x = .data$true_coi, y = .data$estimated_coi)
  ) +
    ggplot2::scale_size_area() +
    ggplot2::geom_abline(color = "red", size = 1) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
      axis.title = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 10)
    ) +
    ggplot2::labs(
      x = "True COI",
      y = "Estimated COI",
      title = sub_title[loop_num]
    )

  # Choose geom_count or geom_boxplot depending on whether we are looking at
  # discrete or continuous data
  if (result_type == "disc") {
    single_plot <- single_plot +
      ggplot2::geom_count(color = "blue", alpha = 0.7, show.legend = FALSE)
  } else if (result_type == "cont") {
    single_plot <- single_plot +
      ggplot2::geom_boxplot(
        color = "blue", alpha = 0.7, show.legend = FALSE,
        ggplot2::aes(group = .data$true_coi)
      )
  }

  single_plot
}

#------------------------------------------------
#' Error plot
#'
#' Creates a plot showing the error of the sensitivity analysis.
#'
#' The function takes in the output of [sensitivity()]. Plots are created using
#' [ggplot2::geom_col()], which creates a simple bar plot. The mean absolute
#' error is plotted in various colors, according to what parameter is being
#' tested. In addition the 95% confidence interval is shown as black vertical
#' lines.
#'
#' @param data The data to be plotted.
#' @param fill The variable the data will be separated by.
#' @param fill_levels The levels for the fill variable.
#' @param title The title of the plot. Default to `NULL`.
#' @param legend_title The text for the legend. Default to `NULL`.
#' @param legend.position The position of the legend. One of `"none"`,
#' `"left"`, `"right"`, `"bottom"`, `"top"`.
#' @param second_fill Indicates if there will be a second fill variable and
#' what it will be.
#'
#' @seealso [ggplot2::geom_col()] for more information on bar plots and the
#' [ggplot2 website](https://ggplot2.tidyverse.org/index.html).
#' @family plotting functions
#' @export
error_plot <- function(data,
                       fill = "coi",
                       fill_levels = NULL,
                       title = NULL,
                       legend_title = fill,
                       legend.position = "right",
                       second_fill = NULL) {

  # Ensure that ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed in order to plot the test.",
      call. = FALSE
    )
  }

  # Check inputs
  assert_in(
    names(data),
    c("predicted_coi", "probability", "param_grid", "boot_error")
  )
  assert_single_string(fill)
  if (!is.null(fill_levels)) assert_string(fill_levels)
  if (!is.null(title)) assert_single_string(title)
  assert_single_string(legend_title)
  assert_single_string(legend.position)
  if (!is.null(second_fill)) assert_single_string(second_fill)

  # Convert data to a tibble
  plot_data <- data$boot_error %>% tidyr::unchop(cols = names(data$boot_error))

  # Make several changes to the data frame. For columns that have more than 1
  # unique value and are not one of c("mae", "lower", "upper", "bias"), convert
  # them to a factor. We also replace 0 with NA for the columns mae, lower, and
  # upper. This makes plotting easier and ensures that error bars are not shown
  # for data points where the mean absolute error is 0.
  plot_data <- plot_data %>%
    dplyr::mutate(dplyr::across(
      where(function(x) dplyr::n_distinct(x) > 1) &
        !tidyselect::all_of(c("mae", "lower", "upper", "bias")),
      as.factor
    )) %>%
    dplyr::mutate(mae = dplyr::na_if(.data$mae, 0)) %>%
    dplyr::mutate(lower = dplyr::na_if(.data$lower, 0)) %>%
    dplyr::mutate(upper = dplyr::na_if(.data$upper, 0))

  # Customize the levels of the fill variable
  if (!is.null(fill_levels)) {
    # Ensure that the number of levels input (fill_levels) are the same as the
    # number of levels for the fill variable
    if (length(fill_levels) != nlevels(plot_data[[fill]])) {
      message <- glue::glue(
        "Number of levels must match:",
        "\n\u2139 Variable has {nlevels(plot_data[[fill]])} levels.",
        "\n\u2716 User specified {length(fill_levels)} levels."
      )
      stop(message, call. = FALSE)
    } else {
      # Customize labels
      levels(plot_data[[fill]]) <- fill_levels
    }
  }

  # Plot the data and return
  error_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$coi, y = .data$mae,
      fill = eval(parse(text = fill))
    )
  ) +
    ggplot2::geom_col(position = "dodge", na.rm = T) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      width = .2,
      position = ggplot2::position_dodge(.9)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = legend.position,
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
      axis.title = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 10)
    ) +
    ggplot2::labs(
      x = "COI", y = "Mean Absolute Error",
      title = title, fill = legend_title
    )

  if (!is.null(second_fill)) {
    error_plot <- error_plot +
      ggplot2::facet_wrap(~ plot_data[[second_fill]], ncol = 2)
  }

  error_plot
}

#------------------------------------------------
#' World map plot
#'
#' Plot a world map showing the COI in each region where reads were sampled
#' from.
#'
#' Creates a world map and overlays the COI in each region. The magnitude of
#' the COI is indicated by both the color and the size of the bubble.
#'
#' @param data The data to be plotted.
#' @param variable The variable the data will plot.
#' @param label The label for the variable.
#' @param alpha The alpha value for the plotted data.
#' @param breaks The breaks for the color scale.
#'
#' @seealso This [website](https://www.r-graph-gallery.com/bubble-map.html) for
#' more information on creating bubble graphs in R.
#' @family plotting functions
#' @export
world_map <- function(data,
                      variable,
                      label = NULL,
                      alpha = 0.1,
                      breaks = c(1, 2)) {

  # Access world map data from ggplot2
  world <- ggplot2::map_data("world")

  # Plot world map
  ggplot2::ggplot() +
    ggplot2::borders("world") +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = .data$long, y = .data$lat, group = .data$group),
      fill = "grey",
      alpha = 0.3
    ) +
    ggplot2::geom_point(
      data = data,
      ggplot2::aes(
        x = .data$long,
        y = .data$lat,
        size = {{ variable }},
        color = {{ variable }}
      ),
      alpha = alpha
    ) +
    ggplot2::scale_colour_viridis_c(
      limits = c(breaks[1], breaks[length(breaks)]),
      breaks = breaks,
      alpha = alpha
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_size(guide = "none") +
    ggplot2::labs(color = label) +
    ggplot2::coord_quickmap(xlim = c(-75, 150), ylim = c(-30, 30))
}
