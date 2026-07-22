library(CrossTalkeR)
library(igraph)
library(colorBlindness)
library(ggplot2)
library(dplyr)

#### From unpublished MF immunity paper

# Plot in- and out-degree barplots for a specific (e.g. filtered/subsetted) igraph CCI network,
# rather than pulling from the pre-computed @rankings table. Useful when the plotted network
# has been modified (edge-weight filtering, added disconnected nodes, etc.) so it no longer
# matches the object's stored rankings.
# weighted = TRUE uses igraph::strength() (sum of edge weights per node, matching CrossTalkeR's
#            own degree ranking behaviour). weighted = FALSE uses igraph::degree() (raw edge counts).
plot_bar_degree_cci <- function(graph, plt_name, path, colors = NULL, weighted = TRUE, top_n = NULL) {
  w <- if (weighted) igraph::E(graph)$weight else NULL

  indegree  <- igraph::strength(graph, mode = "in",  weights = w)
  outdegree <- igraph::strength(graph, mode = "out", weights = w)
  if (!weighted) {
    indegree  <- igraph::degree(graph, mode = "in")
    outdegree <- igraph::degree(graph, mode = "out")
  }

  degree_table <- data.frame(
    nodes = igraph::V(graph)$name,
    InDegree = indegree,
    OutDegree = outdegree
  )

  pdf(paste0(path, plt_name, "_Degree_Barplots.pdf"))
  for (ranking in c("InDegree", "OutDegree")) {
    curr_table <- degree_table %>%
      arrange(get(ranking)) %>%
      as.data.frame()
    if (!is.null(top_n)) {
      curr_table <- tail(curr_table, n = top_n)
    }
    curr_table <- unique(curr_table)
    rownames(curr_table) <- curr_table$nodes

    if (!is.null(colors)) {
      bar_colors <- colors[curr_table$nodes]
      p <- ggplot(curr_table, aes(x = get(ranking), y = reorder(nodes, get(ranking)), fill = nodes)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = bar_colors, guide = "none")
    } else {
      p <- ggplot(curr_table, aes(x = get(ranking), y = reorder(nodes, get(ranking)))) +
        geom_bar(stat = "identity", fill = Blue2DarkOrange18Steps[14])
    }
    print(p +
      ylab("Cell Type") +
      xlab(ranking) +
      ggtitle(paste(plt_name, "-", ranking)) +
      theme_minimal())
  }
  dev.off()

  return(degree_table)
}

# Faithful R port of pyCrossTalkeR's cci_community_layout algorithm.
# The key idea (confirmed from the actual Python source): community detection and layout do NOT
# operate on the raw directed CCI graph. Instead:
#   1) A structural-similarity graph is built: two cell types end up close together if they send/
#      receive interactions to/from similar partners (a symmetrized, degree-normalised similarity
#      kernel over the signed weighted adjacency matrix).
#   2) Communities are detected on that similarity graph (Leiden by default, matching pyCrossTalkeR's
#      use of the leidenalg package with RBConfigurationVertexPartition; Louvain as an alternative).
#   3) Communities are laid out relative to each other (Fruchterman-Reingold, analogous to
#      networkx's spring_layout).
#   4) Nodes are laid out within their own community (same layout algorithm), then compressed
#      towards their community's centre.
# The actual CCI edges (with their real, possibly signed, interaction weights) are drawn separately
# afterwards using these positions - this function only produces the layout, not the plot itself.
# Returns coordinates usable directly as the `coords` argument of new_plot_cci(), plus the detected
# community membership (useful for colouring nodes by community instead of by e.g. lineage).
community_layout_cci <- function(graph,
                                  method = c("leiden", "louvain"),
                                  res = 1.0,
                                  compression = 0.2,
                                  seed = 12,
                                  spring_iterations_community = 1000,
                                  spring_iterations_nodes = 500) {
  method <- match.arg(method)
  set.seed(seed)

  nodes <- igraph::V(graph)$name

  # --- Step 1: structural similarity graph ---
  W <- as.matrix(igraph::as_adjacency_matrix(graph, attr = "weight", sparse = FALSE))
  W <- W[nodes, nodes]

  deg <- rowSums(abs(W))
  deg_inv_sqrt <- ifelse(deg > 0, 1 / sqrt(deg), 0)
  W_deg <- sweep(W, 1, deg_inv_sqrt, `*`)   # scale rows
  W_deg <- sweep(W_deg, 2, deg_inv_sqrt, `*`) # scale columns

  K <- 0.5 * (t(W_deg) %*% W_deg + W_deg %*% t(W_deg))
  diag(K) <- 0
  K[K <= 0] <- 0 # keep only positive similarities, mirroring the Python K_np[i, j] > 0 filter
  dimnames(K) <- list(nodes, nodes)

  g_sim <- igraph::graph_from_adjacency_matrix(K, mode = "upper", weighted = TRUE, diag = FALSE)

  # --- Step 2: community detection (on the similarity graph) ---
  partition <- if (method == "leiden") {
    igraph::cluster_leiden(g_sim,
      objective_function = "modularity",
      resolution = res,
      weights = igraph::E(g_sim)$weight,
      n_iterations = 2
    )
  } else {
    igraph::cluster_louvain(g_sim, weights = igraph::E(g_sim)$weight, resolution = res)
  }
  membership <- igraph::membership(partition)
  communities <- sort(unique(membership))

  # --- Step 3: community-level (macro) layout ---
  contracted <- igraph::contract(g_sim, membership, vertex.attr.comb = "first")
  contracted <- igraph::simplify(contracted, edge.attr.comb = "sum", remove.loops = TRUE)
  meta_layout <- igraph::layout_with_fr(contracted, niter = spring_iterations_community)
  meta_layout <- igraph::norm_coords(meta_layout, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  rownames(meta_layout) <- as.character(communities)

  # --- Step 4: within-community (local) layout, compressed towards the community centre ---
  coords <- matrix(NA_real_, nrow = length(nodes), ncol = 2, dimnames = list(nodes, c("x", "y")))
  for (comm in communities) {
    comm_nodes <- nodes[membership == comm]
    if (length(comm_nodes) == 1) {
      local_layout <- matrix(c(0, 0), nrow = 1)
    } else {
      sub_g <- igraph::induced_subgraph(g_sim, comm_nodes)
      local_layout <- igraph::layout_with_fr(sub_g, niter = spring_iterations_nodes)
      local_layout <- igraph::norm_coords(local_layout, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
    }
    centre <- meta_layout[as.character(comm), ]
    coords[comm_nodes, ] <- sweep(local_layout * compression, 2, centre, FUN = "+")
  }

  return(list(coords = coords, membership = membership))
}

# Helper: assign a distinct colour per detected community, for use as the `colors` argument of
# new_plot_cci() in place of e.g. lineage colours - lets you colour nodes by detected community.
community_colors <- function(membership, palette = NULL) {
  n_comms <- length(unique(membership))
  if (is.null(palette)) {
    palette <- grDevices::rainbow(n_comms)
  }
  node_colors <- palette[membership]
  names(node_colors) <- names(membership)
  return(node_colors)
}


new_plot_cci <- function(graph,
                     colors,
                     plt_name,
                     coords,
                     emax = NULL,
                     leg = FALSE,
                     low = 25,
                     high = 75,
                     ignore_alpha = FALSE,
                     log = FALSE,
                     efactor = 8,
                     vfactor = 12,
                     vnames = T,
                     pg = NULL,
                     vnamescol = NULL, 
                     col_pallet = NULL,
                     standard_node_size = 20,
                     pg_node_size_low = 10,
                     pg_node_size_high = 60,
                     arrow_size = 0.4,
                     arrow_width = 0.8,
                     node_label_position = 1.25,
                     node_label_size = 0.6) {

  # Check Maximal Weight
  if (is.null(emax)) {
    emax <- max(abs(igraph::E(graph)$weight))
  }
  # Using color pallet to up and down regulation
  if(is.null(col_pallet)){
      col_pallet <- colorBlindness::Blue2DarkOrange18Steps
      # Expanding the pallet range
      col_pallet[10] <- '#ffefd7'
  }
  col_pallet <- grDevices::colorRampPalette(col_pallet)(201)
  # Checking looops
  edge_start <- igraph::ends(graph,
                             es = igraph::E(graph),
                             names = FALSE)
  # Scale nodes coordinates
  if (nrow(coords) != 1) {
    coords_scale <- scale(coords)
  } else {
    coords_scale <- coords
  }
  # It will make the loops in a correct angle
loop_angle <- atan2(
  coords_scale[igraph::V(graph)$name, 2],
  coords_scale[igraph::V(graph)$name, 1]
)
  # Setting node colors
  igraph::V(graph)$color <- colors[igraph::V(graph)$name]
  ## Color scheme
  we <- round(
    oce::rescale(igraph::E(graph)$weight,
                 xlow = (-emax),
                 xhigh = emax,
                 rlow = 1,
                 rhigh = 200,
                 clip = TRUE),
    0)
  igraph::E(graph)$color <- col_pallet[we]
  alpha_cond <- (igraph::E(graph)$inter > low) & (igraph::E(graph)$inter < high)
  alpha <- ifelse(alpha_cond, 0, igraph::E(graph)$inter)
  subgraph <- igraph::delete.edges(graph,
                                   igraph::E(graph)[alpha == 0 | is.na(alpha)]
  )
  if (!ignore_alpha) {
    igraph::E(graph)$color <- scales::alpha(igraph::E(graph)$color, alpha)
  }
  ## Thickness and arrow size
  if (is.null(pg)) {
    igraph::V(graph)$size <- standard_node_size
  } else {
    igraph::V(graph)$size <- oce::rescale(pg,
                                          xlow = quantile(x = pg, prob = 0.25),
                                          xhigh = quantile(x = pg, prob = 0.75),
                                          rlow = pg_node_size_low,
                                          rhigh = pg_node_size_high,
                                          clip = TRUE)
  }

  if (log) {
    igraph::E(graph)$width <- ifelse(igraph::E(graph)$inter != 0,
                                     log2(1 + igraph::E(graph)$inter),
                                     0) * efactor
  } else {
    igraph::E(graph)$width <- ifelse(igraph::E(graph)$inter != 0,
                                     igraph::E(graph)$inter,
                                     0) * efactor
  }
  igraph::E(graph)$arrow.size <- arrow_size
  igraph::E(graph)$arrow.width <- igraph::E(graph)$width + arrow_width
  igraph::E(graph)$loop.angle <- NA
  if (sum(edge_start[, 2] == edge_start[, 1]) != 0) {
    igraph::E(graph)$loop.angle[which(edge_start[, 2] == edge_start[, 1])] <- loop_angle[edge_start[which(edge_start[, 2] == edge_start[, 1]), 1]]
    igraph::E(graph)$loop.angle[which(edge_start[, 2] != edge_start[, 1])] <- 0
  }
  coords_scale[, 1] <- scales::rescale(coords_scale[, 1], from = c(-1, 1), to = c(-2, 2))
  coords_scale[, 2] <- scales::rescale(coords_scale[, 2], from = c(-1, 1), to = c(-2, 2))
  plot(graph,
       layout = coords_scale,
       xlim = c(-4, 4),
       ylim = c(-4, 4),
       rescale = F,
       edge.curved = 0.5,
       vertex.label = NA,
       vertex.shape = "circle",
       margin = 0.0,
       loop.angle = igraph::E(graph)$loop.angle,
       edge.label = NA,
       main = plt_name
  )

  x <- coords_scale[, 1] * node_label_position
  y <- coords_scale[, 2] * node_label_position
  coord_ratio <- coords_scale[, 1] / coords_scale[, 2]
  angle <- ifelse(
    atan(-coord_ratio) * (180 / pi) < 0,
    90 + atan(-coord_ratio) * (180 / pi),
    270 + atan(-coord_ratio) * (180 / pi))
  if (vnames) {
    if (!is.null(vnamescol)) {
      for (i in seq_len(length(x))) {
        graphics::text(x = x[i],
                      y = y[i],
                      labels = igraph::V(graph)$name[i],
                      adj = NULL,
                      pos = NULL,
                      cex = 0.6,
                      col = vnamescol[igraph::V(graph)$name[i]],
                      xpd = TRUE)
      }
    } else {
      for (i in seq_len(length(x))) {
        graphics::text(x = x[i],
                      y = y[i],
                      labels = igraph::V(graph)$name[i],
                      adj = NULL,
                      pos = NULL,
                      cex = 0.6,
                      col = "black",
                      xpd = TRUE)
      }

    }
  }

  if (leg) {
    # Thicknesse legend
    amin <- min(igraph::E(graph)$inter[igraph::E(graph)$inter != 0])
    amax <- max(igraph::E(graph)$inter)
    e_wid_sp <- c(amin,
                  amin + amax / 2,
                  amax)
    graphics::legend("topleft",
                    legend = round(e_wid_sp, 1),
                    col = "black",
                    title = "Percentage of the interactions",
                    pch = NA,
                    bty = "n",
                    cex = 1,
                    lwd = e_wid_sp,
                    lty = c(1, 1, 1),
                    horiz = FALSE)

    v <- igraph::V(graph)$size
    # Pagerank legend
    if (!is.null(pg)) {
      a <- graphics::legend('bottomleft',
                            title = "Node Pagerank",
                            legend = c("", "", ""),
                            pt.cex = c(min(v) + 1, mean(v), max(v)) / 12, col = 'black',
                            pch = 21, pt.bg = 'black', box.lwd = 0, y.intersp = 2)
      graphics::text(a$rect$left + a$rect$w, a$text$y,
                    c(round(min(pg), 2), round(mean(pg), 2), round(max(pg), 2)), pos = 2)
    }
    # Edge Colormap
    if (min(igraph::E(graph)$weight) < 0 & max(igraph::E(graph)$weight) > 0) {
      leg <- netdiffuseR::drawColorKey(seq(1, 200),
                                       tick.marks = c(1, 101, 200),
                                       color.palette = col_pallet,
                                       labels = c(-round(emax, 3), 0, round(emax, 3)),
                                       nlevels = 200,
                                       main = "Weights",
                                       pos = 2,
                                       key.pos = c(0.98, 1.0, 0.0, 0.2),
                                       border = "transparent")
    } else if (max(igraph::E(graph)$weight) < 0) {
      leg <- netdiffuseR::drawColorKey(seq(1, 100),
                                       tick.marks = c(1, 101),
                                       color.palette = col_pallet[1:101],
                                       labels = c(-round(emax, 3), 0),
                                       nlevels = 100,
                                       main = "Weights",
                                       pos = 2,
                                       key.pos = c(0.98, 1.0, 0.0, 0.2),
                                       border = "transparent")
    }else {
      leg <- netdiffuseR::drawColorKey(seq(100, 200),
                                       tick.marks = c(100, 200),
                                       color.palette = col_pallet[100:201],
                                       labels = c(0, round(emax, 3)),
                                       nlevels = 100,
                                       main = "Weights",
                                       pos = 2,
                                       key.pos = c(0.98, 1.0, 0.0, 0.2),
                                       border = "transparent")
    }
  }
}
