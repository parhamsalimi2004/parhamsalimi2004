import requests
import pandas as pd

CBIOPORTAL_BASE_URL = "https://www.cbioportal.org/api"

def load_tcga_mutations(study_id = "acc_tcga"):
    molecular_profile_id = f"{study_id}_mutations"
    sample_list_id = f"{study_id}_all"

    url = f"{CBIOPORTAL_BASE_URL}/molecular-profiles/{molecular_profile_id}/mutations/fetch"
    params = {"projection": "DETAILED"}
    body = {"sampleListId": sample_list_id}

    response = requests.post(url, params=params, json=body, timeout=60)
    response.raise_for_status()
    mutations = response.json()

    df = pd.json_normalize(mutations)
    return df

def load_cosmic_labels(cmc_path, passenger_sample_size=50000, random_state=0):
    cmc = pd.read_csv(cmc_path, sep="\t", compression="gzip", low_memory=False)
# only mutations with significance tier 1 & 2
    drivers = cmc[
        cmc["MUTATION_SIGNIFICANCE_TIER"].isin(["1","2"])
    ].copy()
    drivers["label"]=1 # label driver mutations
# only mutations with significance tier "Other"
# only mutations with NaN CGC tier
    passengers = cmc[
        (cmc["MUTATION_SIGNIFICANCE_TIER"] == "Other") &
        (cmc["CGC_TIER"].isna())
    ].copy()
    passengers["label"] = 0 # label passenger(negative) mutations

    if passenger_sample_size is not None and len(passengers) > passenger_sample_size:
        passengers = passengers.sample(
            n=passenger_sample_size,
            random_state=random_state
        )
    labels = pd.concat([drivers, passengers], ignore_index=True) #stack objects

    keep_columns = [
        "GENE_NAME",
        "Mutation AA",
        "MUTATION_SIGNIFICANCE_TIER",
        "CGC_TIER",
        "COSMIC_SAMPLE_TESTED",
        "COSMIC_SAMPLE_MUTATED",
        "label",
    ]
    return labels[keep_columns]

def load_string_network(links_path, info_path, min_confidence=400):
    info = pd.read_csv(info_path, sep="\t")
    id_to_gene = dict(zip(info["#string_protein_id"], info["preferred_name"]))
    links = pd.read_csv(links_path, sep=" ")
    links = links[links["combined_score"]>= min_confidence]
    edges = pd.DataFrame({
        "geneA": links["protein1"].map(id_to_gene),
        "geneB": links["protein2"].map(id_to_gene)
    })
    edges = edges.dropna(subset=["geneA","geneB"])

    return edges

if __name__ == "__main__":
    df = load_tcga_mutations()
    print(f"pulled {len(df)} mutations across {df['patientId'].nunique() if 'patientId' in df.columns else '?'} patients")
    print(df.head())