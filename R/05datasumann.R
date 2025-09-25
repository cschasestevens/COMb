#' Annotation Summary
#'
#' Summarizes annotation information for a specific dataset.
#'
#' @param ld A data list containing annotations assigned by ms_input().
#' @param md A vector of annotation levels to display in the plot.
#' @return List containing a formatted annotation table and sunburst plot
#' summarizing annotation information.
#' @examples
#'
#' # ms_plot_anno(
#' #   d,
#' #   c("Major.Class", "Class", "Subclass")
#' # )
#'
#' @export
ms_plot_anno <- function(ld, md) {
  # Load dataset
  d <- ld
  md1 <- md
  # Extract annotations
  d_anno <- d[["data"]][["anno"]][, md1]
  # Create Nodes for selected levels
  d_anno <- as.data.frame(
    lapply(
      seq.int(1, ncol(d_anno), 1),
      function(x) {
        factor(d_anno[[x]], levels = sort(unique(d_anno[[x]])))
      }
    )
  )
  ## Node counts for each level
  ct1 <- dplyr::bind_rows(
    lapply(
      seq.int(1, ncol(d_anno), 1),
      function(x) {
        c1 <- dplyr::count(
          d_anno,
          d_anno[[x]]
        )
        names(c1) <- c("Node", "N")
        return(c1)
      }
    )
  )
  ct1 <- setNames(
    aggregate(
      ct1[["N"]],
      list(ct1[["Node"]]),
      function(x) sum(x)
    ),
    c("Node", "N")
  )
  ## Node proportions
  ct2 <- dplyr::bind_rows(
    lapply(
      seq.int(0, ncol(d_anno) - 1, 1),
      function(x) {
        c1 <- dplyr::mutate(
          dplyr::count(
            d_anno,
            d_anno[, 1:(x + 1)]
          ),
          Proportion = round((n / sum(n)) * 100, digits = 2), # nolint
        )
        c1 <- setNames(
          c1[, c((ncol(c1) - 2), ncol(c1))],
          c("Node", "Proportion")
        )
        return(c1)
      }
    )
  )
  ct2 <- setNames(
    aggregate(
      ct2[["Proportion"]],
      list(ct2[["Node"]]),
      function(x) sum(x)
    ),
    c("Node", "Proportion")
  )
  ## Level
  ct3 <- dplyr::bind_rows(
    lapply(
      seq.int(1, ncol(d_anno), 1),
      function(x) {
        c1 <- as.data.frame(levels(d_anno[[x]]))
        c1[["Level"]] <- x
        names(c1) <- c("Node", "Level")
        return(c1)
      }
    )
  )
  ct3 <- ct3[!duplicated(ct3), ]
  ## Class Label
  ct4 <- dplyr::bind_rows(
    lapply(
      seq.int(0, ncol(d_anno) - 1, 1),
      function(x) {
        c1 <- dplyr::mutate(
          dplyr::count(
            d_anno,
            d_anno[, 1:(x + 1)]
          ),
          Proportion = round((n / sum(n)) * 100, digits = 2), # nolint
        )
        c1 <- setNames(
          c1[, c((ncol(c1) - 2), 1)],
          c("Node", "Class")
        )
        return(c1)
      }
    )
  )
  ct4 <- ct4[!duplicated(ct4[["Node"]]), ]
  node_anno <- dplyr::left_join(
    ct1,
    dplyr::left_join(
      ct2,
      dplyr::left_join(
        ct3,
        ct4,
        by = "Node"
      ),
      by = "Node"
    ),
    by = "Node"
  )
  ## Node label positions
  an1 <- unlist(
    lapply(
      seq.int(1, ncol(d_anno), 1),
      function(x) unique(as.character(d_anno[[x]]))
    )
  )
  node_anno[["Node"]] <- factor(node_anno[["Node"]], levels = an1)
  node_anno <- node_anno[
    order(
      node_anno[["Class"]],
      node_anno[["Node"]]
    ),
  ]
  node_anno <- dplyr::mutate(
    dplyr::group_by(
      node_anno,
      .data[["Level"]] # nolint
    ),
    End = 2 * pi * cumsum(N) / sum(N), # nolint
    Start = dplyr::lag(End, default = 0), # nolint
    Mid = 0.5 * (Start + End) # nolint
  )
  ## Formatted node labels
  node_anno[["Label"]] <- paste(
    node_anno[["Node"]], ", ",
    node_anno[["N"]], ", ",
    node_anno[["Proportion"]], "%",
    sep = ""
  )
  node_anno[["Label"]] <- ifelse(
    node_anno[["Proportion"]] < quantile(node_anno[["Proportion"]])[[4]],
    "",
    node_anno[["Label"]]
  )
  # Edges
  ed1 <- dplyr::bind_rows(
    lapply(
      seq.int(1, ncol(d_anno) - 1, 1),
      function(x) {
        d1 <- setNames(
          data.frame(
            d_anno[[x]],
            d_anno[[x + 1]]
          ),
          c("from", "to")
        )
        return(d1)
      }
    )
  )
  ed1 <- ed1[!duplicated(ed1), ]

  # Final Network
  net1 <- igraph::graph_from_data_frame(
    ed1,
    vertices = node_anno
  )

  # Sunburst Plot
  pnet1 <- ggraph::ggraph(
    net1,
    layout = "partition",
    circular = TRUE
  ) +
    ggraph::geom_node_arc_bar(
      ggplot2::aes(
        fill = .data[["Class"]], # nolint
        r0 = ((Level - 1) / ncol(d_anno)) ^ 1.25, # nolint
        r = ((Level) / ncol(d_anno)) ^ 1.25,
        start = Start, # nolint
        end = End, # nolint
        alpha = Proportion # nolint
      ),
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(values = col_univ()) + # nolint
    shadowtext::geom_shadowtext(
      ggplot2::aes(
        x = (((Level) / ncol(d_anno)) ^ 1.25 +
            ((Level - 1) / ncol(d_anno)) ^ 1.25
        ) * sin(Mid) / 2, # nolint
        y = (((Level) / ncol(d_anno)) ^ 1.25 +
            ((Level - 1) / ncol(d_anno)) ^ 1.25
        ) * cos(Mid) / 2,
        label = Label # nolint
      ),
      bg.color = "white",
      color = "black",
      alpha = 0.8,
      bg.r = 0.2
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      legend.key.size = ggplot2::unit(1, "cm"),
      legend.title = ggplot2::element_text(
        size = 18,
        face = "bold",
        vjust = 0.2
      ),
      plot.margin = ggplot2::unit(c(1.5, 1.5, 1.5, 1.5), "cm"),
      legend.position = c("right"),
      legend.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(
        size = 20,
        face = "bold",
        hjust = 0.5
      )
    )
  return(
    list(
      "table" = node_anno,
      "plot" = pnet1
    )
  )
}
