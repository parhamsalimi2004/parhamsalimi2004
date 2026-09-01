import pandas as pd
import numpy as np
import networkx as nx

def compute_hotspot_recurrence(labels):
    result = labels.copy()
    result["recurrence_count"] = result["COSMIC_SAMPLE_MUTATED"].fillna(0)
    result["log_recurrence_count"] = np.log1p(result["recurrence_count"])
    return result

def assign_recurrence_group(features, n_recurrent_groups=2):
    result = features.copy()
    is_singleton = result["recurrence_count"] <= 1
    result["recurrence_group"] = "singleton" # initially label everything as singleton
    recurrent_vals = result.loc[~is_singleton, "recurrence_count"]

    if len(recurrent_vals) > 0:
        ranked = recurrent_vals.rank(method="first")
        n_groups = min(n_recurrent_groups, len(recurrent_vals))
        sub_groups = pd.qcut(ranked, q=n_groups, labels=False, duplicates="drop")
        result.loc[~is_singleton, "recurrence_group"] = "recurrent_" + sub_groups.astype(str)
    
    return result
    
def compute_network_centrality(gene_list, string_edges):
    G = nx.from_pandas_edgelist(string_edges, "geneA", "geneB")
    degree = nx.degree_centrality(G)
    betweenness = nx.betweenness_centrality(G)

    return pd.DataFrame({
        "gene": gene_list,
        "degree_centrality": [degree.get(g,0.0) for g in gene_list],
        "betweenness_centrality": [betweenness.get(g,0.0) for g in gene_list]
    })

def compute_cancer_gene_proximity(gene_list: list, string_edges: pd.DataFrame, known_cancer_genes: set, smoothing_strength=10):
    # fraction of neighbour cancer genes
    # distance to nearest cancer gene
    G = nx.from_pandas_edgelist(string_edges, "geneA", "geneB")
    known_in_graph = set(known_cancer_genes) & set(G.nodes())
    global_rate = len(known_in_graph) / G.number_of_nodes()

    G_temp = G.copy()
    G_temp.add_edges_from([("__SUPER_SOURCE__", g) for g in known_in_graph])
    raw_distances = nx.single_source_shortest_path_length(G_temp, "__SUPER_SOURCE__")

    fractions, neighbor_counts, distances = [], [], []
    for gene in gene_list:
        if gene not in G:
            fractions.append(global_rate)  # no evidence at all -> just use the average
            neighbor_counts.append(0)
            distances.append(np.nan)
            continue
        neighbors = set(G.neighbors(gene))
        n_total = len(neighbors)
        n_cancer = len(neighbors & known_in_graph)
        smoothed = (n_cancer + smoothing_strength * global_rate) / (n_total + smoothing_strength)
        fractions.append(smoothed)
        neighbor_counts.append(n_total)
        distances.append(raw_distances.get(gene, np.nan) - 1 if gene in raw_distances else np.nan)

    return pd.DataFrame({
        "gene": gene_list,
        "fraction_neighbors_cancer_genes": fractions,
        "neighbor_count": neighbor_counts,
        "distance_to_nearest_cancer_gene": distances,
    })
