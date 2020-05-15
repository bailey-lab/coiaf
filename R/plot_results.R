#------------------------------------------------
#' @title Function to plot results from \code{\link{coi_test}}
#'
#' @description Plots a result
#'
#' @param data The data to be plotted
#' @param plot_dims A list representing the number of rows and columns
#' @param change_param The parameter being changed
#' @param change_param_val The values the changed parameter takes
#'
#' @return Plot
#'
#' @export

sensitivity_plot <- function(data,
                        plot_dims,
                        change_param = NULL,
                        change_param_val = NULL){

  plot_df = data$predicted_coi %>%
    tidyr::gather(true_COI, estimated_COI) %>%
    tidyr::extract(true_COI, c("true_COI", "loop_number"), "coi_(.+)_(.+)") %>%
    dplyr::mutate_all(as.numeric)

  num_loops <- unique(plot_df$loop_number)
  myplots <- lapply(num_loops,
                    sensitivity_plot_element,
                    data = plot_df,
                    change_param = change_param,
                    change_param_val = change_param_val)

  final_plot <- ggpubr::ggarrange(plotlist = myplots,
                                  labels = "AUTO",
                                  nrow = plot_dims[1],
                                  ncol = plot_dims[2])

  return(final_plot)
}

#------------------------------------------------
#' @title Helper function for \code{\link{sensitivity_plot}}
#'
#' @description Plots a singular result
#'
#' @param data The data to be plotted
#' @param loop_num The loop number
#' @inheritParams sensitivity_plot
#'
#' @return Plot
#'
#' @keywords internal

sensitivity_plot_element <- function(data, loop_num, change_param, change_param_val){

  ggplot2::ggplot(dplyr::filter(data, loop_number == loop_num),
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
}
