#------------------------------------------------
#' @title Plot sensitivity analysis
#'
#' @description The function is used to plot any sensitivity analysis that is
#' run. It takes in the output of \code{\link{coi_test}} and creates a grid of
#' plots.
#'
#' @param data The data to be plotted.
#' @param plot_dims A list representing the number of rows and columns our plots
#' will be split into.
#' @param change_param The title of the plot. The title specifies the parameter
#' that is being changed.
#' @param change_param_val The values the changed parameter ranges over.
#' @param title The title of the figure
#' @param caption The caption of the figure
#'
#' @return A grid of plots that represents the sensitivity analysis
#'
#' @export

sensitivity_plot <- function(data,
                             plot_dims,
                             change_param = NULL,
                             change_param_val = NULL,
                             title = NULL,
                             caption = NULL){

  # Ensure that ggplot2 and ggpubr are installed
  if (!requireNamespace("ggplot2", quietly = TRUE) &
      !requireNamespace("ggpubr", quietly = TRUE)) {
    stop('Packages \"ggplot2\" and \"ggpubr\" must be installed in order to plot the tests.',
         call. = FALSE)
  } else if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop('Package \"ggplot2\" must be installed in order to plot the tests.',
         call. = FALSE)
  } else if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop('Package \"ggpubr\" must be installed in order to plot the tests.',
         call. = FALSE)
  }

  # Check inputs
  assert_eq(names(data), c("predicted_coi", "param_grid", "error_bias"))
  assert_pos_int(plot_dims, zero_allowed = FALSE)
  assert_length(plot_dims, 2)
  if (!is.null(change_param)) {assert_string(change_param)}
  if (!is.null(change_param_val)) {assert_vector(change_param_val)}
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
  suggested_dims = plot_dims[1] * plot_dims[2]
  needed_dims    = length(num_loops)
  if (!all(suggested_dims >= needed_dims)) {
    stop(sprintf("Not enough panels have been specified. User input %s panels, but %s panels needed.",
                    suggested_dims, needed_dims), call. = FALSE)
  }

  # We then call a helper function: sensitivity_plot_element, that creates each
  # individual plot and store these plots as a list
  myplots <- lapply(num_loops,
                    sensitivity_plot_element,
                    data = plot_df,
                    change_param = change_param,
                    change_param_val = change_param_val)

  # Arrange the plots
  if (plot_dims[1] * plot_dims[2] == 1){
    # Do not include a panel label if there is only one plot
    arranged_plots <- ggpubr::ggarrange(plotlist = myplots,
                                        nrow = plot_dims[1],
                                        ncol = plot_dims[2])
  } else {
    # Include panel labels if there are more than one plots
    arranged_plots <- ggpubr::ggarrange(plotlist = myplots,
                                        labels = "AUTO",
                                        font.label = list(size = 10),
                                        nrow = plot_dims[1],
                                        ncol = plot_dims[2])
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
#' @return A single panel of the sensitivity analysis.
#'
#' @keywords internal

sensitivity_plot_element <- function(data,
                                     loop_num,
                                     change_param,
                                     change_param_val){

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
    ggplot2::labs(x = "True COI",
                  y = "Estimated COI",
                  title = paste0(change_param, change_param_val[loop_num]))

  return(single_plot)
}

#------------------------------------------------
#' @title Plot the mean absolute errors for a test
#'
#' @description Plots the mean absolute error and the confidence interval.
#'
#' @param data The data to be plotted.
#' @param fill The variable the data will be seperated by.
#' @param fill_levels The levels for the fill variable.
#' @param title The title of the plot. Default to \code{NULL}.
#' @param legend_title The text for the legend. Default to \code{NULL}.
#' @param legend.position The position of the legend. One of \code{"none",
#' "left", "right", "bottom", "top"}.
#'
#' @return The plot
#'
#' @export

error_plot <- function(data, fill = "COI", fill_levels = NULL, title = NULL,
                       legend_title = fill, legend.position = "bottom"){

  # Ensure that ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop('Package \"ggplot2\" must be installed in order to plot the test.',
         call. = FALSE)
  }

  # Check inputs
  assert_single_string(fill)
  if (!is.null(fill_levels)) {assert_string(fill_levels)}
  if (!is.null(title)) {assert_single_string(title)}
  assert_single_string(legend_title)
  assert_single_string(legend.position)

  # Convert data to a tibble
  plot_data <- data$boot_error %>% tidyr::unnest(cols = names(data$boot_error))

  # Make several changes to the data frame. For columns that have more than 1
  # unique value and are not one of c("mae", "lower", "upper", "bias"), convert
  # them to a factor. We also replace 0 with NA for the columns mae, lower, and
  # upper. This makes plotting easier and ensures that error bars are not shown
  # for data points where the mean absolute error is 0.
  plot_data <- plot_data %>%
    dplyr::mutate(dplyr::across(where(function(x) dplyr::n_distinct(x) > 1) &
                                  !all_of(c("mae", "lower", "upper", "bias")),
                                as.factor)) %>%
    dplyr::mutate(mae   = dplyr::na_if(mae, 0)) %>%
    dplyr::mutate(lower = dplyr::na_if(lower, 0)) %>%
    dplyr::mutate(upper = dplyr::na_if(upper, 0))

  # Customize the levels of the fill variable
  if (!is.null(fill_levels)){
    # Ensure that the number of levels input (fill_levels) are the same as the
    # number of levels for the fill variable
    if (length(fill_levels) != nlevels(plot_data[[fill]])){
      stop(strwrap(sprintf("The length of the fill_levels variable (%s) does not
      match with the number of levels of the fill variable (%s)",
                           length(fill_levels), nlevels(plot_data[[fill]]))),
           call. = FALSE)
    }
    else {
      # Customize labels
      levels(plot_data[[fill]]) <- fill_levels
    }
  }

  # Plot the data and return
  error_plot <- ggplot2::ggplot(plot_data,
                                ggplot2::aes(x = COI, y = mae,
                                             fill = eval(parse(text = fill)))) +
    ggplot2::geom_col(position = position_dodge(), na.rm = T) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                           width = .2, position = position_dodge(.9)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = legend.position,
                   plot.title   = ggplot2::element_text(hjust = 0.5, size = 10),
                   axis.title   = ggplot2::element_text(size = 7),
                   legend.title = ggplot2::element_text(size = 7),
                   legend.text  = ggplot2::element_text(size = 7)) +
    ggplot2::labs(y = "Mean Absolute Error", title = title, fill = legend_title)

  return(error_plot)

}
