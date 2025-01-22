plot_pangenome_cumulative = function (fit, plot = TRUE, legend = TRUE, text_size = 14, color_pallete = 6, 
          facet = FALSE, smooth = FALSE) 
{
  if (class(fit) != "panfit") {
    purrr::map(fit, ~{
      if (class(.x) != "panfit") 
        stop("fit is not of class `panfit`!")
      validate_panfit(.x)
    })
  }
  else {
    validate_panfit(fit)
    fit <- list(fit)
  }
  plot_data <- purrr::imap_dfr(fit, ~{
    temp_tree <- .x$tree
    temp_tree$edge.length <- .x$data$acc
    return(tibble::tibble(pangenome = .y, acc = ape::node.depth.edgelength(temp_tree), 
                          core = ape::node.depth.edgelength(.x$tree), branch = ifelse(c(1:(.x$tree$Nnode + 
                                                                                             length(.x$tree$tip.label))) < length(.x$tree$tip.label), 
                                                                                      "terminal", "internal")))
  })
  if (!plot) {
    return(plot_data)
  }
  if (length(fit) > 1) {
    gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$core, 
                                                  y = .data$acc, colour = .data$pangenome)) + ggplot2::geom_point(ggplot2::aes(shape = .data$branch))
    if (facet) {
      plot_data$pangenome2 = plot_data$pangenome
      gg = ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$core, y = .data$acc, colour = .data$pangenome)) + 
        ggplot2::geom_point(data = plot_data[,2:5], ggplot2::aes(x = .data$core, y = .data$acc), colour = "grey", size = 0.5) +
        ggplot2::geom_point()
      gg <- gg + ggplot2::facet_wrap(~pangenome, ncol = round(sqrt(length(unique(plot_data$pangenome)))))
    }
  }
  else {
    gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$core, 
                                                  y = .data$acc)) + ggplot2::geom_point(ggplot2::aes(shape = .data$branch))
  }
  if (smooth) {
    warning("Adding trend line using 'ggplot2::geom_smooth'. This is not the panstripe fit!")
    gg <- gg + ggplot2::geom_smooth(method = "glm", method.args = list(family = "quasipoisson"), 
                                    level = 0.95)
  }
  gg <- gg + ggplot2::scale_colour_brewer(type = "qual", palette = color_pallete) + 
  ggplot2::scale_fill_brewer(type = "qual", palette = color_pallete) +
  ggplot2::theme_bw(base_size = text_size) + ggplot2::xlab("cumulative core branch distance") + 
  ggplot2::ylab("cumulative genomic divergence events")
  if (!legend) {
    gg <- gg + ggplot2::theme(legend.position = "none") +
      ggplot2::guides(color = "none", fill = "none")
  }
  gg
  return(gg)
}
