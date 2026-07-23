library(CrossTalkeR)
library(igraph)
library(colorBlindness)
library(ggplot2)
library(dplyr)
library(ggraph)
source("./module_03_CrossTalker_extra_plotting_new.R")

object_Jak2 <- readRDS("../data/MF/MF_JAK2_LR_data_final.Rds")
subgraph_MF_pre_Jak2 <- subgraph.edges(object_Jak2@graphs$MF_pre_treatment, 
                        E(object_Jak2@graphs$MF_pre_treatment)[E(object_Jak2@graphs$MF_pre_treatment)$weight > 40])
nodes_Jak2 <- V(subgraph_MF_pre_Jak2)$name

# Lineage Colors
meta_data_table <- read.csv("../data/MF/MF_subcluster_subset_metadata.csv", row.names = 1)
meta_data_table <- meta_data_table %>%
  select(subcluster, lineage_broad)
test_table <- data.frame(lineage = meta_data_table$lineage_broad, cluster = meta_data_table$subcluster)
test_table <- test_table %>%
  unique() %>%
  group_by(lineage) %>%
  summarize(new_col <- list(cluster))

niche_HSPCs <- sort(append(unlist(as.vector(test_table[[1, 2]])), paste(unlist(as.vector(test_table[[1, 2]])), "Mut", sep = "_")))
niche_Lymphoid <- sort(append(unlist(as.vector(test_table[[2, 2]])), paste(unlist(as.vector(test_table[[2, 2]])), "Mut", sep = "_")))
niche_Lymphoid_pro <- sort(append(unlist(as.vector(test_table[[3, 2]])), paste(unlist(as.vector(test_table[[3, 2]])), "Mut", sep = "_")))
niche_Myeloid <- sort(append(unlist(as.vector(test_table[[4, 2]])), paste(unlist(as.vector(test_table[[4, 2]])), "Mut", sep = "_")))
niche_Myeloid_pro <- sort(append(unlist(as.vector(test_table[[5, 2]])), paste(unlist(as.vector(test_table[[5, 2]])), "Mut", sep = "_")))
niche_Non_heamatopoietic <- sort(append(unlist(as.vector(test_table[[6, 2]])), paste(unlist(as.vector(test_table[[6, 2]])), "Mut", sep = "_")))

nichescol <- c(
  rep("#dec234", length(niche_HSPCs)),
  rep("#006639", length(niche_Lymphoid)),
  rep("#ff648e", length(niche_Myeloid_pro)),
  rep("#aa0051", length(niche_Myeloid)),
  rep("#b1db8e", length(niche_Lymphoid_pro)),
  rep("#019ce5", length(niche_Non_heamatopoietic))
)

niches <- c(niche_HSPCs,
            niche_Lymphoid,
            niche_Myeloid,
            niche_Myeloid_pro,
            niche_Lymphoid_pro,
            niche_Non_heamatopoietic)
names(nichescol) <- niches

all_nodes<-nodes_Jak2

add_not_connected_nodes <- function(subgraph, node_list){
  nodes_filtered <- V(subgraph)$name
  new_vertices <- setdiff(node_list, nodes_filtered)
  subgraph <- add_vertices(subgraph, length(new_vertices), name = new_vertices)
  return(subgraph)
}

template <- igraph::make_full_graph(
  n = length(all_nodes),
  directed = TRUE,
  loops = TRUE
)
smaller_circle <- igraph::layout.circle(template)

fil_nodes <- niches[niches %in% all_nodes]
rownames(smaller_circle) <- fil_nodes

# Plot
subgraph_MF_pre_Jak2 <- add_not_connected_nodes(subgraph_MF_pre_Jak2, all_nodes)
pdf("../figures/revision/MF/CCI_MF.pdf", width = 4.5, height = 4.5)
new_plot_cci(subgraph_MF_pre_Jak2,
                  paste0("MF pre-treatment mutational interactions", " CCI"),
                  leg = FALSE,
                  low = 0,
                  high = 0 / 100,
                  ignore_alpha = FALSE,
                  log = TRUE,
                  efactor = 2,
                  vfactor = 12,
                  vnames = TRUE,
                  pg = object_Jak2@rankings[["MF_pre_treatment"]]$Pagerank[V(subgraph_MF_pre_Jak2)$name], 
                  vnamescol = NULL, 
                  colors = nichescol[V(subgraph_MF_pre_Jak2)$name],
                  coords = smaller_circle[V(subgraph_MF_pre_Jak2)$name, ],,
                  col_pallet = NULL, 
                  standard_node_size = 20,
                  pg_node_size_low = 10, 
                  pg_node_size_high = 40,
                  arrow_size = 0.3, 
                  arrow_width = 0.6, 
                  node_label_position = 1.25,
                  node_label_size = 0.6)
dev.off()

degree_Jak2 <- plot_bar_degree_cci(
  graph = subgraph_MF_pre_Jak2,
  plt_name = "MF_pre_Jak2",
  path = "../figures/revision/MF/",
  colors = nichescol,
  weighted = TRUE   # set FALSE for raw connection counts instead of weighted strength
)

## Get CCI table for network
edge_list <- igraph::as_data_frame(subgraph_MF_pre_Jak2, what = "edges")
CCI_table <- edge_list[, c("from", "to", "weight")]

write.csv(
    CCI_table,
    "../res_data/revision/MF/Jak2_CCI_table.csv",
    row.names = FALSE
)

comm_0_2=c('SPP1_Macrophages', 'MK_prog', 'SPP1_Macrophages_Mut', 'CD16_monocytes_Mut', 'Classical_monocytes_Mut', 'Inflammatory_CD14_monocytes_Mut', 'OLCs', 'Adipo_CAR')
pdf('../figures/revision/F7/F7F.pdf', 18, 5)


plot_sankey(lrobj_tbl = object_Jak2@tables$MF_pre_treatment, threshold=25, ligand_cluster=comm_0_2, receptor_cluster=comm_0_2)


dev.off()