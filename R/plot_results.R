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
#'
#' @return A grid of plots that represents the sensitivity analysis
#'
#' @export

sensitivity_plot <- function(data,
                             plot_dims,
                             change_param = NULL,
                             change_param_val = NULL){
  # Check inputs
  assert_eq(names(data), c("predicted_coi", "param_grid", "error_bias"))
  assert_pos_int(plot_dims, zero_allowed = FALSE)
  assert_length(plot_dims, 2)
  if (!assert_null(change_param)) {assert_string(change_param)}
  if (!assert_null(change_param_val)) {assert_vector(change_param_val)}

  # Convert the predicted_coi dataframe into a long format. More specifically,
  # establish a column for the true coi and a column for the estimated coi
  # value. In some cases, we change parameters other than the COI. When this is
  # the case, we want to make multiple plots and see the effect of the changing
  # parameter. Thus, we establish another column: loop number, that tells us
  # how often the other parameter is changing. In the end, there will be three
  # columns: true_COI, estimated_COI, and loop_number.
  plot_df = data$predicted_coi %>%
    tidyr::gather(true_COI, estimated_COI) %>%
    tidyr::extract(true_COI, c("true_COI", "loop_number"), "coi_(.+)_(.+)") %>%
    dplyr::mutate_all(as.numeric)

  # We determine how many different panels there will be by finding the unique
  # values of loop_number
  num_loops <- unique(plot_df$loop_number)

  # Need to make sure that there are enough panels to include all the graphs
  assert_greq(plot_dims[1] * plot_dims[2], length(num_loops))

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
                                        nrow = plot_dims[1],
                                        ncol = plot_dims[2])
  }

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

  single_plot <- ggplot2::ggplot(dplyr::filter(data, loop_number == loop_num),
                                 aes(x = true_COI, y = estimated_COI)) +
    geom_count(color = "blue", alpha = 0.7, show.legend = FALSE) +
    geom_abline(color = "red") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)) +
    labs(x = "True COI",
         y = "Estimated COI",
         title = paste(change_param, change_param_val[loop_num], sep = " = "))

  return(single_plot)
}
