import pandas as pd

cmc = pd.read_csv("data/CancerMutationCensus_AllData_Tsv_v104_GRCh37/CancerMutationCensus_AllData_v104_GRCh37.tsv.gz", sep="\t", compression="gzip", low_memory=False)
cgc = pd.read_csv("data/example_grch37/Cosmic_CancerGeneCensus_v104_GRCh37.tsv", sep="\t")

#print(cmc.columns.tolist())
#print(cgc.columns.tolist())

#print(cmc['MUTATION_SIGNIFICANCE_TIER'].value_counts(dropna=False))
#print(cmc['CGC_TIER'].value_counts(dropna=False))

#from data_sources import load_cosmic_labels
#labels = load_cosmic_labels("data/CancerMutationCensus_AllData_Tsv_v104_GRCh37/CancerMutationCensus_AllData_v104_GRCh37.tsv.gz")
#print(result.shape)
#print(result['label'].value_counts())

#print(labels[labels["label"] == 1]["COSMIC_SAMPLE_MUTATED"].describe())
#print(labels["COSMIC_SAMPLE_MUTATED"].describe())

#from features import compute_hotspot_recurrence, assign_recurrence_group

#feats = compute_hotspot_recurrence(labels)
#grouped = assign_recurrence_group(feats)
#print(grouped["recurrence_group"].value_counts())
#print(grouped.groupby("recurrence_group")["label"].mean())

#print(cmc["MUTATION_SIGNIFICANCE_TIER"].value_counts(dropna=False))
#print(cmc["MUTATION_SIGNIFICANCE_TIER"].dtype)


from data_sources import load_cosmic_labels
#from features import compute_hotspot_recurrence, assign_recurrence_group

labels = load_cosmic_labels("data/CancerMutationCensus_AllData_Tsv_v104_GRCh37/CancerMutationCensus_AllData_v104_GRCh37.tsv.gz")

#feats = compute_hotspot_recurrence(labels)
#print(feats[["recurrence_count", "log_recurrence_count"]].describe())

#grouped = assign_recurrence_group(feats)
#print(grouped["recurrence_group"].value_counts())
#print(grouped.groupby("recurrence_group")["label"].mean())

#protinfo = pd.read_csv("data/9606.protein.info.v12.0.txt",sep='\t')
#protlink = pd.read_csv("data/9606.protein.links.v12.0.txt",sep=' ')

#print(protinfo.columns.to_list())
#print(protinfo.head())
#print(protlink.columns.to_list())
#print(protlink.head())

from data_sources import load_string_network
#from features import compute_network_centrality
from features import compute_cancer_gene_proximity
edges = load_string_network("data/9606.protein.links.v12.0.txt", "data/9606.protein.info.v12.0.txt")
#print(edges.shape)
known_cancer_genes = set(cgc["GENE_SYMBOL"])

gene_list = labels["GENE_NAME"].unique().tolist()
proximity = compute_cancer_gene_proximity(gene_list, edges, known_cancer_genes)
print(proximity.sort_values("fraction_neighbors_cancer_genes", ascending=False).head(10))
#centrality = compute_network_centrality(gene_list, edges)
#print(centrality.sort_values("degree_centrality", ascending=False).head(10))

