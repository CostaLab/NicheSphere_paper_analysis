import pandas as pd
import pycrosstalker
import numpy as np 
import networkx as nx
import igraph as ig
import leidenalg 
from matplotlib.colors import TwoSlopeNorm
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.patches import Patch
from community_layout.layout_class import CommunityLayout
import nichesphere
import scanpy as sc

### from https://github.com/CostaLab/pyCrossTalkeR/blob/dev/pycrosstalker/plots/plot.py
def get_community_cmap(n):
    return plt.get_cmap("tab20", n)

def cci_community_layout(
    CCI_table,
    method="leiden",
    res=1.0,
    compression=0.2,
    seed=12,
    spring_k=0.1,
    spring_iterations_community=1000,
    spring_iterations_nodes=500,
    edge_width_scale=3.0,
    edge_alpha_intra=0.7,
    edge_alpha_inter=0.7,
    node_size=300,
    arrow_size=20,
    figsize=[10, 10],
    return_nw=False
):
    """
    Unified community layout for cell-cell communication networks.
    Parameters
    ----------
    method : str
        Community detection algorithm: "leiden" or "louvain".
    spring_iterations_nodes : int
        Only used when method="leiden" (per-community spring layout iterations).
    spring_k : float
        Spring layout k parameter. Leiden default: 0.1, Louvain default: 75.
    """
    # ---- Build similarity matrix ----
    all_cells = pd.Index(sorted(set(CCI_table["from"]) | set(CCI_table["to"])))
    W = CCI_table.pivot_table(
        index="from", columns="to", values="weight", fill_value=0.0
    ).reindex(index=all_cells, columns=all_cells, fill_value=0.0)
    cells = W.index.to_numpy()
    #W = CCI_table.pivot_table(
    #    index="from", columns="to", values="weight", fill_value=0.0
    #)
    #cells = W.index.to_numpy()

    deg = W.abs().sum(axis=1)
    deg_inv_sqrt = 1.0 / np.sqrt(deg.replace(0, np.nan))
    W_deg = W.mul(deg_inv_sqrt, axis=0).mul(deg_inv_sqrt, axis=1).fillna(0)
    K = 0.5 * (W_deg.T @ W_deg + W_deg @ W_deg.T)
    K_np = K.values

    # ---- Build undirected similarity graph ----
    G = nx.Graph()
    G.add_nodes_from(cells)
    for i in range(len(cells)):
        for j in np.argsort(K_np[i])[::-1]:
            if j != i and K_np[i, j] > 0:
                G.add_edge(cells[i], cells[j], weight=K_np[i, j])

    # ---- Community detection + layout ----
    if method == "leiden":
        node_to_idx = {n: i for i, n in enumerate(G.nodes())}
        idx_to_node = {i: n for n, i in node_to_idx.items()}
        edges = [
            (node_to_idx[u], node_to_idx[v], d["weight"])
            for u, v, d in G.edges(data=True)
        ]

        g = ig.Graph(
            n=len(G.nodes()),
            edges=[(u, v) for u, v, _ in edges],
            edge_attrs={"weight": [w for _, _, w in edges]},
            directed=False,
        )
        #g = ig.Graph(
        #    edges=[(u, v) for u, v, _ in edges],
        #    edge_attrs={"weight": [w for _, _, w in edges]},
        #    directed=False,
        #)
        partition = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            weights="weight",
            resolution_parameter=res,
            seed=seed,
        )
        for idx, comm in enumerate(partition.membership):
            G.nodes[idx_to_node[idx]]["community"] = comm

        communities = [
            {n for n, d in G.nodes(data=True) if d["community"] == cid}
            for cid in sorted(set(partition.membership))
        ]

        # Community-level graph for macro layout
        G_comm = nx.Graph()
        for cid in range(len(communities)):
            G_comm.add_node(cid)
        for u, v, d in G.edges(data=True):
            cu, cv = G.nodes[u]["community"], G.nodes[v]["community"]
            if cu != cv:
                w = d.get("weight", 1.0)
                if G_comm.has_edge(cu, cv):
                    G_comm[cu][cv]["weight"] += w
                else:
                    G_comm.add_edge(cu, cv, weight=w)

        pos_comm = nx.spring_layout(
            G_comm, k=spring_k, iterations=spring_iterations_community, seed=seed
        )
        pos = {}
        for cid, nodes in enumerate(communities):
            subG = G.subgraph(nodes)
            local_pos = nx.spring_layout(
                subG, k=spring_k, iterations=spring_iterations_nodes, seed=seed
            )
            center = pos_comm[cid]
            for n, p in local_pos.items():
                pos[n] = center + compression * p

    elif method == "louvain":
        cl = CommunityLayout(
            G,
            community_compression=compression,
            layout_algorithm=nx.spring_layout,
            layout_kwargs={"k": spring_k, "iterations": spring_iterations_community},
            community_algorithm=nx.algorithms.community.louvain_communities,
            community_kwargs={"resolution": res, "seed": seed, "weight": "weight"},
        )
        pos = cl.full_positions
        communities = list(cl.communities())

    else:
        raise ValueError(f"Unknown method '{method}'. Choose 'leiden' or 'louvain'.")

    ###
    missing = set(G.nodes()) - set(pos.keys())
    if missing:
        rng = np.random.default_rng(seed)
        center = np.mean(list(pos.values()), axis=0) if pos else np.zeros(2)
        for n in missing:
            pos[n] = center + rng.normal(scale=0.05, size=2)
    ###

    # ---- Build directed signed graph for drawing ----
    G_signed_dir = nx.DiGraph()
    for _, row in CCI_table.iterrows():
        if method == "louvain" and row["weight"] == 0:
            continue
        G_signed_dir.add_edge(row["from"], row["to"], weight=row["weight"])

    weights = np.array([d["weight"] for _, _, d in G_signed_dir.edges(data=True)])
    norm = TwoSlopeNorm(
        vmin=-abs(weights).max(), vcenter=0.0, vmax=abs(weights).max()
    )
    edge_colors = cm.coolwarm(norm(weights))
    edge_widths = 0.6 + edge_width_scale * np.abs(weights) / np.max(np.abs(weights))

    node_comm = {n: cid for cid, nodes in enumerate(communities) for n in nodes}
    n_comms = len(communities)
    cmap = get_community_cmap(n_comms)
    nodes = list(G_signed_dir.nodes())
    node_colors = [cmap(node_comm[n]) for n in nodes]
    edge_alphas = [
        edge_alpha_intra if node_comm[u] == node_comm[v] else edge_alpha_inter
        for u, v in G_signed_dir.edges()
    ]

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=figsize)

    nx.draw_networkx_nodes(
        G_signed_dir, pos, nodelist=nodes, node_color=node_colors,
        node_size=node_size, edgecolors="black", linewidths=0.7, ax=ax,
    )
    for (u, v), c, w, a in zip(G_signed_dir.edges(), edge_colors, edge_widths, edge_alphas):
        nx.draw_networkx_edges(
            G_signed_dir, pos, edgelist=[(u, v)], edge_color=[c], width=w,
            alpha=a, arrows=True, arrowstyle="-|>", arrowsize=arrow_size,
            connectionstyle="arc3,rad=0.15", ax=ax,
        )

    texts = [
        ax.text(
            x, y, node, fontsize=9, ha="center", va="center",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.5, pad=0.2),
            zorder=10,
        )
        for node, (x, y) in pos.items()
    ]
    adjust_text(
        texts, ax=ax, expand_points=(1.2, 1.2), expand_text=(1.2, 1.2),
        force_text=0.8, force_points=0.2,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.6),
    )

    comm_handles = [
        Patch(facecolor=cmap(i), edgecolor="black", label=f"Community {i}")
        for i in range(n_comms)
    ]
    ax.legend(
        handles=comm_handles, title="Communication communities",
        bbox_to_anchor=(1.02, 1.0), loc="upper left",
    )

    sm = cm.ScalarMappable(cmap=cm.coolwarm, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.4)
    cbar.set_label("Differential LR interaction score")

    ax.set_title(
        f"Cell–cell communication community network\n"
        f"{method.capitalize()} resolution: {res}",
        fontsize=14,
    )
    ax.axis("off")
    plt.tight_layout()
    plt.show()
    if return_nw:
        cols = cmap(np.linspace(0, 1, 4))
        clist=[mcolors.to_hex(color) for color in cols]
        nichesphere.tl.assign_properties(G_signed_dir, 
                                 communities=cl.communities(), 
                                 colors=clist, 
                                 pos=pos, 
                                 simmilarity_weights=False, 
                                 node_coord_sf=1000,
                                 #g_unsigned=G
                                 )
        return G_signed_dir
    else:
        return fig


def plot_degree_barplots(CCI_table, weight_col="weight", figsize=(14, 6), sort_by="total_n"):
    """
    Barplots of in- and out-degree (edge counts) and in- and out-strength
    (summed weights) per cell type in a CCI table.

    Parameters
    ----------
    CCI_table : pd.DataFrame
        Must contain columns 'from', 'to', and `weight_col`.
    weight_col : str
        Name of the column holding interaction weights/scores.
    figsize : tuple
        Figure size passed to plt.subplots.
    sort_by : str
        One of {"total_n", "total_weight", "name"} - how to order the bars.

    Returns
    -------
    fig : matplotlib.figure.Figure
    deg : pd.DataFrame
        Per-cell-type table of out_n, out_weight, in_n, in_weight, total_n, total_weight.
    """
    out_deg = (CCI_table.groupby("from")[weight_col]
               .agg(["count", "sum"])
               .rename(columns={"count": "out_n", "sum": "out_weight"}))
    in_deg = (CCI_table.groupby("to")[weight_col]
              .agg(["count", "sum"])
              .rename(columns={"count": "in_n", "sum": "in_weight"}))

    # outer join so cells that only send or only receive still show up (fill 0)
    deg = out_deg.join(in_deg, how="outer").fillna(0.0)
    deg["total_n"] = deg["out_n"] + deg["in_n"]
    deg["total_weight"] = deg["out_weight"] + deg["in_weight"]

    if sort_by == "name":
        deg = deg.sort_index()
    else:
        deg = deg.sort_values(sort_by, ascending=False)

    x = np.arange(len(deg))
    width = 0.35

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    axes[0].bar(x - width / 2, deg["out_n"], width, label="Out-degree (sending)")
    axes[0].bar(x + width / 2, deg["in_n"], width, label="In-degree (receiving)")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(deg.index, rotation=90)
    axes[0].set_ylabel("Number of interactions")
    axes[0].set_title("Interaction count (in vs out)")
    axes[0].legend()

    axes[1].bar(x - width / 2, deg["out_weight"], width, label="Out-strength (sending)")
    axes[1].bar(x + width / 2, deg["in_weight"], width, label="In-strength (receiving)")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(deg.index, rotation=90)
    axes[1].set_ylabel("Sum of interaction weight")
    axes[1].set_title("Interaction weight (in vs out)")
    axes[1].legend()

    plt.tight_layout()
    return fig, deg


def _build_digraph(CCI_table, weight_col="weight"):
    """
    Build a directed graph from a CCI table the same way cci_community_layout does,
    so degree centrality is computed on the same network you're visualizing elsewhere.
    """
    G = nx.DiGraph()
    for _, row in CCI_table.iterrows():
        G.add_edge(row["from"], row["to"], weight=row[weight_col])
    return G


def plot_degree_centrality(CCI_table, weight_col="weight", figsize=(7, 6), sort=True):
    """
    Horizontal barplots of in-degree and out-degree centrality (networkx),
    as two SEPARATE figures.

    Uses nx.in_degree_centrality / nx.out_degree_centrality, which are UNWEIGHTED:
    each is (number of distinct neighbors) / (n_nodes - 1). Interaction weights
    are ignored. If you want interaction weight to factor in, use a weighted
    variant instead (e.g. dict(G.in_degree(weight="weight"))).

    Parameters
    ----------
    CCI_table : pd.DataFrame with columns 'from', 'to', weight_col
    sort : bool - sort bars by centrality value (descending, so highest is at top)

    Returns
    -------
    (fig_in, fig_out) : matplotlib Figures
    (in_cent, out_cent) : dicts of {node: centrality}
    """
    G = _build_digraph(CCI_table, weight_col)

    in_cent = nx.in_degree_centrality(G)
    out_cent = nx.out_degree_centrality(G)

    def _make_plot(cent_dict, title, color):
        s = pd.Series(cent_dict)
        if sort:
            s = s.sort_values(ascending=True)  # ascending so highest bar ends up on top after barh
        else:
            s = s.sort_index(ascending=False)

        fig, ax = plt.subplots(figsize=figsize)
        ax.barh(s.index, s.values, color=color)
        ax.set_xlabel("Degree centrality")
        ax.set_title(title)
        plt.tight_layout()
        #ax.spines.set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #plt.axis('off') 
        return fig

    fig_in = _make_plot(in_cent, "In-degree centrality", "tab:orange")
    fig_out = _make_plot(out_cent, "Out-degree centrality", "tab:blue")

    return (fig_in, fig_out), (in_cent, out_cent)


def _barh_plot(values_dict, title, xlabel, color, figsize=(7, 6), sort=True):
    s = pd.Series(values_dict)
    if sort:
        s = s.sort_values(ascending=True)  # ascending so highest ends up on top after barh
    else:
        s = s.sort_index(ascending=False)

    fig, ax = plt.subplots(figsize=figsize)
    ax.barh(s.index, s.values, color=color)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    plt.tight_layout()
    return fig


def plot_betweenness_centrality(CCI_table, weight_col="weight", use_weight=False,
                                 figsize=(7, 6), sort=True):
    """
    Horizontal barplot of betweenness centrality (nx.betweenness_centrality).

    use_weight : bool
        If True, edge 'weight' is used as a DISTANCE in shortest-path computations
        (networkx convention: higher weight = "farther"/less likely to be on a
        shortest path). Since your interaction weights mean the OPPOSITE
        (higher = stronger interaction), this is almost never what you want as-is.
        Left False by default (topology only). If you do want weight to count,
        consider passing a distance column equal to 1/weight instead.

    Returns
    -------
    fig : matplotlib Figure
    values : dict of {node: betweenness centrality}
    """
    G = _build_digraph(CCI_table, weight_col)
    values = nx.betweenness_centrality(G, weight=weight_col if use_weight else None)
    fig = _barh_plot(values, "Betweenness centrality", "Betweenness centrality",
                      "tab:green", figsize, sort)
    return fig, values


def plot_pagerank(CCI_table, weight_col="weight", use_weight=True,
                   figsize=(7, 6), sort=True):
    """
    Horizontal barplot of PageRank (nx.pagerank).

    use_weight : bool
        If True (default), edge 'weight' is used to weight transition
        probability along edges (higher weight = more "flow" toward that
        neighbor) - this matches the usual meaning of interaction strength,
        unlike betweenness's distance convention above.

    Returns
    -------
    fig : matplotlib Figure
    values : dict of {node: pagerank score}
    """
    G = _build_digraph(CCI_table, weight_col)
    values = nx.pagerank(G, weight=weight_col if use_weight else None)
    fig = _barh_plot(values, "PageRank", "PageRank score", "tab:purple", figsize, sort)
    return fig, values