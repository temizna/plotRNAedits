#' Plot 192 trinucleotide RNA edit profile
#'
#' Plot relative contribution of 192 trinucleotides
#' This function is a variation of plot_96_proile from MutationalPatterns
#' @param mut_matrix 192 trinucleotide profile matrix
#' @param ymax Y axis maximum value, default = 0.2
#' @param colors Optional 12 value color vector.
#' @param condensed More condensed plotting format. Default = F.
#' @return 192 trinucleotide profile plot
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#'
#' @examples
#' ## See the 'edit_matrix()' example for how we obtained the
#' ## RNA edit matrix information:
#' edit_mat <- load(system.file("data/edit_mat.rda"))
#'
#' ## Plot the 192-profile of two samples
#' plot_192_edit_profile(edit_mat, ymax=0.15)
#'
#'
#' @export

plot_192_edit_profile<-function (mut_matrix, colors = NA, ymax = 0.2, condensed = FALSE)
{
  freq <- full_context <- substitution <- context <- NULL
  if (is.na(colors)) {
    colors <- c("yellow","black","red","cyan","black","red",
                "cyan","yellow","red","cyan","yellow","black")
  }
  if (length(colors) != 12) {
    stop("Provide colors vector with length 6", call. = FALSE)
  }
  norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
  tb <- norm_mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>%
    dplyr::mutate(substitution = stringr::str_replace(full_context,
                                                      "\\w\\[(.*)\\]\\w", "\\1"), context = stringr::str_replace(full_context,
                                                                                                                 "\\[.*\\]", "\\.")) %>% dplyr::select(-full_context) %>%
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample",
                        values_to = "freq") %>% dplyr::mutate(sample = factor(sample,
                                                                              levels = unique(sample)))
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.6
    spacing <- 0.5
  }
  plot <- ggplot(data = tb, aes(x = context, y = freq, fill = substitution,
                                width = width)) + geom_bar(stat = "identity", colour = "black",
                                                           linewidth = 0.2) + scale_fill_manual(values = colors) + facet_grid(sample ~
                                                                                                                           substitution) + ylab("Relative contribution") + coord_cartesian(ylim = c(0,
                                                                                                                                                                                                    ymax)) + scale_y_continuous(breaks = seq(0, ymax, 0.1)) +
    guides(fill = "none") + theme_bw() + theme(axis.title.y = element_text(size = 10,
                                                                           vjust = 1), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 8),
                                               axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5),
                                               strip.text.x = element_text(size = 8), strip.text.y = element_text(size = 8),
                                               panel.grid.major.x = element_blank(), panel.spacing.x = unit(spacing,
                                                                                                            "lines"))
  return(plot)
}
